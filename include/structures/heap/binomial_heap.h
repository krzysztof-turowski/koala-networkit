#ifndef KOALA_BINOMIALHEAP_H
#define KOALA_BINOMIALHEAP_H

#include <cassert>
#include <functional>
#include <iterator>

namespace Koala {
	template <class Key, class Compare = std::less<Key> >
	class BinomialHeap {
		struct node {
			node *parent, *child, *next;
			unsigned degree;
			Key key;

			inline node(const Key &key = Key()) : parent(0), child(0), next(0), degree(0), key(key) { }

			inline void insert(node *A) {
				A->parent = this, A->next = child, child = A, degree++;
			}

			void clear() {
				if (next) {
					next->clear();
        }
				if (child) {
					child->clear();
        }
				delete this;
			}
		};

		node *root, *minimum;
		unsigned nodes;
		Compare function;
	public:
		typedef Key value_type;

		class iterator : public std::iterator<std::input_iterator_tag, node> {
			friend class BinomialHeap<Key, Compare>;

			node *data;
		public:
			typedef typename BinomialHeap<Key, Compare>::value_type value_type;
			typedef value_type* pointer;
			typedef value_type& reference;

			iterator(node *data = 0) : data(data) { }
			iterator(const iterator &A) : data(A.data) { }
			bool operator==(const iterator &A)	{ return data == A.data; }
			bool operator!=(const iterator &A)	{ return data != A.data; }
			reference operator*() const			{ return data->key; }
			pointer operator->() const			{ return &**this; }
		};

		inline BinomialHeap(const Compare& = Compare());
		~BinomialHeap();

		const value_type& top() const;
		iterator push(const value_type&);
		void pop();

		void update(iterator, const value_type&);
		void erase(iterator);

		void clear();

		unsigned size() const;
		bool empty() const;
	private:
		inline node* join(node*, node*);
		inline node* reverse(node*);
		inline node* cut(node*);
	};

	template <class Key, class Compare>
	inline BinomialHeap<Key, Compare>::BinomialHeap(const Compare &function) : root(0), minimum(0), nodes(0), function(function) { }

	template <class Key, class Compare>
	BinomialHeap<Key, Compare>::~BinomialHeap()
	{
		clear();
	}

	template <class Key, class Compare>
	const typename BinomialHeap<Key, Compare>::value_type& BinomialHeap<Key, Compare>::top() const
	{
		return minimum->key;
	}

	template <class Key, class Compare>
	typename BinomialHeap<Key, Compare>::iterator BinomialHeap<Key, Compare>::push(const typename BinomialHeap<Key, Compare>::value_type &key)
	{
		nodes++;
		node *A = new node(key);

		if(root == 0)
			return iterator(root = minimum = A);

		root = join(root, A);
		if(!function(A->key, minimum->key))
			minimum = A;

		while(minimum->parent)
			minimum = minimum->parent;
		return iterator(A);
	}

	template <class Key, class Compare>
	void BinomialHeap<Key, Compare>::pop()
	{
		nodes--;

		if(root == minimum)
			root = root->next;
		else
		{
			node *A = root;
			while(A->next != minimum)
				A = A->next;
			A->next = minimum->next;
		}

		if(nodes == 0)
		{
			delete minimum, minimum = 0;
			return;
		}

		node *child = minimum->child;
		if(child)
		{
			for(node *A = child; A; A = A->next)
				A->parent = 0;
			root = root ? join(root, reverse(child)) : reverse(child);
		}

		delete minimum, minimum = root;
		if(minimum)
			for(node *A = root->next; A; A = A->next)
				if(!function(A->key, minimum->key))
					minimum = A;
	}

	template <class Key, class Compare>
	void BinomialHeap<Key, Compare>::update(typename BinomialHeap<Key, Compare>::iterator element, const typename BinomialHeap<Key, Compare>::value_type &key)
	{
		node *A = element.data;
		assert(!function(key, A->key));

		A->key = key;
		if(!function(key, minimum->key))
			minimum = A;

		if(!A->parent || !function(A->parent->key, A->key))
		{
			while(minimum->parent)
				minimum = minimum->parent;
			return;
		}

		node *start = 0, *previous = 0, *B = A, *C = A->parent, *D;
		while(C)
		{
			D = C->child, C->child = B->next;
			if(B == D)
				D = C, C->degree--;
			else
			{
				node *E = D;
				while(B != E->next)
					E = E->next, C->degree--;
				E->next = C, C->degree -= 2;
			}
			B->next = start, B = C, C = C->parent, start = previous, previous = D;
		}

		if(B == root)
			root = root->next;
		else
		{
			C = root;
			while(B != C->next)
				C = C->next;
			C->next = B->next;
		}
		B->next = start, start = previous;

		if(start)
		{
			for(B = start; B; B = B->next)
				B->parent = 0;
			root = root ? join(root, reverse(start)) : reverse(start);
		}
		A->parent = 0, A->next = 0, root = root ? join(root, A) : A;

		while(minimum->parent)
			minimum = minimum->parent;
	}

	template <class Key, class Compare>
	void BinomialHeap<Key, Compare>::erase(iterator element)
	{
		nodes--;

		node *A = element.data;
		if(nodes == 0)
		{
			delete A, root = minimum = 0;
			return;
		}

		node *start = A->child, *previous = A->child, *next, *B = A, *C = A->parent;
		while(C)
		{
			next = C->child, C->child = B->next;
			if(B == next)
				next = C, C->degree--;
			else
			{
				node *D = next;
				while(B != D->next)
					D = D->next, C->degree--;
				D->next = C, C->degree -= 2;
			}
			B->next = start, B = C, C = C->parent, start = previous, previous = next;
		}

		if(B == root)
			root = root->next;
		else
		{
			C = root;
			while(B != C->next)
				C = C->next;
			C->next = B->next;
		}
		B->next = start, start = previous;

		if(start)
		{
			for(B = start; B; B = B->next)
				B->parent = 0;
			root = root ? join(root, reverse(start)) : reverse(start);
		}

		if(minimum == A)
		{
			minimum = root;
			for(B = root->next; B; B = B->next)
				if(!function(B->key, minimum->key))
					minimum = B;
		}
		delete A;
	}

	template <class Key, class Compare>
	void BinomialHeap<Key, Compare>::clear()
	{
		if(root)
			root->clear(), root = 0, nodes = 0;
	}

	template <class Key, class Compare>
	unsigned BinomialHeap<Key, Compare>::size() const
	{
		return nodes;
	}

	template <class Key, class Compare>
	bool BinomialHeap<Key, Compare>::empty() const
	{
		return root == 0;
	}

	template <class Key, class Compare>
	inline typename BinomialHeap<Key, Compare>::node* BinomialHeap<Key, Compare>::join(node *A, node *B)
	{
		node *start, *C;
		if(A->degree <= B->degree)
			start = C = A, A = A->next;
		else
			start = C = B, B = B->next;
		while(A && B)
		{
			if(A->degree <= B->degree)
				C->next = A, A = A->next;
			else
				C->next = B, B = B->next;
			C = C->next;
		}
		C->next = A ? A : B;

		for(A = 0, B = start, C = B->next; C; C = B->next)
			if(B->degree != C->degree || (C->next && C->degree == C->next->degree))
				A = B, B = C;
			else if(!function(B->key, C->key))
				B->next = C->next, B->insert(C);
			else
			{
				if(A)
					A->next = C;
				else
					start = C;
				C->insert(B), B = C;
			}
		return start;
	}

	template <class Key, class Compare>
	inline typename BinomialHeap<Key, Compare>::node* BinomialHeap<Key, Compare>::reverse(node *A)
	{
		node *B = A->next, *C;
		A->next = 0;
		while(B)
			C = B->next, B->next = A, A = B, B = C;
		return A;
	}

	template <class Key, class Compare>
	inline typename BinomialHeap<Key, Compare>::node* BinomialHeap<Key, Compare>::cut(node *A)
	{
		node *B = A->next, *C;
		A->next = 0;
		while(B)
			C = B->next, B->next = A, A = B, B = C;
		return A;
	}
}

#endif
