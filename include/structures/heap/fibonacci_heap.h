#ifndef KOALA_FIBONACCIHEAP_H
#define KOALA_FIBONACCIHEAP_H

#include <cassert>
#include <cstring>
#include <functional>
#include <iterator>

namespace Koala
{
	template <class Key, class Compare = std::less<Key> >
	class FibonacciHeap
	{
		struct node
		{
			node *parent, *child, *previous, *next;
			unsigned flag;
			Key key;

			inline node(const Key &key = Key()) : parent(0), child(0), flag(0), key(key)
			{
				previous = next = this;
			}

			inline void insert(node *A)
			{
				next->previous = A->previous, A->previous->next = next, next = A, A->previous = this;
			}

			inline void remove()
			{
				previous->next = next, next->previous = previous, previous = next = this;
			}

			void clear()
			{
				if(child)
					child->clear();
				if(previous != this)
					next->previous = previous, next->clear();
				delete this;
			}
		};

		node *root, **degrees;
		unsigned nodes;
		Compare function;
	public:
		typedef Key value_type;

		class iterator : public std::iterator<std::input_iterator_tag, node>
		{
			friend class FibonacciHeap<Key, Compare>;

			node *data;
		public:
			typedef typename FibonacciHeap<Key, Compare>::value_type value_type;
			typedef value_type* pointer;
			typedef value_type& reference;

			iterator(node *data = 0) : data(data) { }
			iterator(const iterator &A) : data(A.data) { }
			bool operator==(const iterator &A)	{ return data == A.data; }
			bool operator!=(const iterator &A)	{ return data != A.data; }
			reference operator*() const			{ return data->key; }
			pointer operator->() const			{ return &**this; }
		};

		inline FibonacciHeap(const Compare& = Compare());
		~FibonacciHeap();

		const value_type& top() const;
		iterator push(const value_type&);
		void pop();

		void update(iterator, const value_type&);
		void erase(iterator);
		void clear();

		unsigned size() const;
		bool empty() const;
	};

	template <class Key, class Compare>
	inline FibonacciHeap<Key, Compare>::FibonacciHeap(const Compare &function) : root(0), nodes(0), function(function)
	{
		unsigned size = sizeof(unsigned) << 3;
		degrees = new node*[size], memset(degrees, 0, size * sizeof(node*));
	}

	template <class Key, class Compare>
	FibonacciHeap<Key, Compare>::~FibonacciHeap()
	{
		clear(), delete[] degrees;
	}

	template <class Key, class Compare>
	const typename FibonacciHeap<Key, Compare>::value_type& FibonacciHeap<Key, Compare>::top() const
	{
		return root->key;
	}

	template <class Key, class Compare>
	typename FibonacciHeap<Key, Compare>::iterator FibonacciHeap<Key, Compare>::push(const typename FibonacciHeap<Key, Compare>::value_type &key)
	{
		nodes++;

		node *A = new node(key);
		if(!root)
			return iterator(root = A);

		root->insert(A);
		if(!function(A->key, root->key))
			root = A;
		return iterator(A);
	}

	template <class Key, class Compare>
	void FibonacciHeap<Key, Compare>::pop()
	{
		nodes--;

		node *A = root->child, *B;
		if(A)
		{
			B = A;
			do
			{
				B->parent = 0, B = B->next;
			}
			while(A != B);
			root->insert(A);
		}

		if(nodes == 0)
		{
			root = 0, delete root;
			return;
		}

		node **degrees = this->degrees, *C;
		unsigned degree_max = 0, degree;
		for(A = root->next, B = A->next; A != root; degrees[degree] = A, A = B, B = A->next)
		{
			while(degrees[degree = A->flag >> 1])
			{
				C = degrees[degree];
				if(!function(C->key, A->key))
					C = A, A = degrees[degree];
				degrees[degree] = 0, C->remove(), C->parent = A, C->flag &= ~1;
				if(A->child)
					A->flag += 2, A->child->insert(C);
				else
					A->flag = 2, A->child = C;
			}

			if(degree > degree_max)
				degree_max = degree;
		}
		root->remove(), delete root;

		for(degree = 0; degree <= degree_max; degree++)
			if(degrees[degree])
			{
				root = degrees[degree], degrees[degree] = 0, degree++;
				break;
			}
		for(; degree <= degree_max; degree++)
			if(degrees[degree])
			{
				if(!function(degrees[degree]->key, root->key))
					root = degrees[degree];
				degrees[degree] = 0;
			}
	}

	template <class Key, class Compare>
	void FibonacciHeap<Key, Compare>::update(const typename FibonacciHeap<Key, Compare>::iterator element, const typename FibonacciHeap<Key, Compare>::value_type &key)
	{
		node *A = element.data;
		assert(!function(key, A->key));

		A->key = key;
		node *B = A->parent;
		if(!B)
		{
			if(!function(key, root->key))
				root = A;
			return;
		}
		else if(function(key, B->key))
			return;

		while(1)
		{
			if(A == A->next)
				B->child = 0;
			else
			{
				if(A == B->child)
					B->child = A->next;
				A->remove(), A->flag &= ~1;
			}
			B->flag -= 2, root->insert(A), A->parent = 0;
			if(!function(A->key, root->key))
				root = A;

			if(!B->parent)
				return;
			if(!(B->flag & 1))
			{
				B->flag |= 1;
				return;
			}
			A = B, B = B->parent;
		}
	}

	template <class Key, class Compare>
	void FibonacciHeap<Key, Compare>::erase(const typename FibonacciHeap<Key, Compare>::iterator element)
	{
		node *A = element.data, *B = A->parent, *C = A;
		if(!B)
		{
			root = A, pop();
			return;
		}

		while(1)
		{
			if(A == A->next)
				B->child = 0;
			else
			{
				if(A == B->child)
					B->child = A->next;
				A->remove(), A->flag &= ~1;
			}
			B->flag -= 2, root->insert(A), A->parent = 0;

			if(!B->parent)
				break;
			if(!(B->flag & 1))
			{
				B->flag |= 1;
				break;
			}
			A = B, B = B->parent;
		}

		root = C, pop();
	}

	template <class Key, class Compare>
	void FibonacciHeap<Key, Compare>::clear()
	{
		if(root)
			root->clear(), root = 0, nodes = 0;
	}

	template <class Key, class Compare>
	unsigned FibonacciHeap<Key, Compare>::size() const
	{
		return nodes;
	}

	template <class Key, class Compare>
	bool FibonacciHeap<Key, Compare>::empty() const
	{
		return !root;
	}
}

#endif
