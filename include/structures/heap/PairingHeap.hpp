#ifndef KOALA_PAIRINGHEAP_H
#define KOALA_PAIRINGHEAP_H

#include <cassert>
#include <functional>
#include <iterator>

namespace Koala
{
	template <class Key, class Compare = std::greater<Key> >
	class PairingHeap
	{
		class node
		{
			friend class PairingHeap<Key, Compare>;

			node *parent, *child, *previous, *next;
			unsigned degree;
		public:
			Key key;
		private:
			inline node(const Key&key = Key()) : parent(0), child(0), previous(0), next(0), degree(0), key(key) { }

			inline void insert(node *A)
			{
				if(child)
					child->previous = A;
				A->parent = this, A->previous = 0, A->next = child, child = A, degree++;
			}

			inline void remove()
			{
				if(this == parent->child)
					parent->child = next;
				else
					previous->next = next;
				if(next)
					next->previous = previous;
				parent->degree--, parent = previous = next = 0;
			}

			void clear()
			{
				if(next)
					next->clear();
				if(child)
					child->clear();

				delete this;
			}
		};

		node *root;
		unsigned nodes;
		Compare function;
	public:
		typedef Key value_type;

		class iterator : public std::iterator<std::input_iterator_tag, node>
		{
			friend class PairingHeap<Key, Compare>;

			node *data;
		public:
			typedef typename PairingHeap<Key, Compare>::value_type value_type;
			typedef value_type* pointer;
			typedef value_type& reference;

			iterator(node *data = 0) : data(data) { }
			iterator(const iterator &A) : data(A.data) { }
			bool operator==(const iterator &A)	{ return data == A.data; }
			bool operator!=(const iterator &A)	{ return data != A.data; }
			reference operator*() const			{ return data->key; }
			pointer operator->() const			{ return &**this; }
		};

		inline PairingHeap(const Compare& = Compare());
		~PairingHeap();

		value_type& top() const;
		iterator push(const value_type&);
		void pop();

		void update(iterator, const value_type&);
		void erase(iterator);

		void clear();

		unsigned size() const;
		bool empty() const;
	};

	template <class Key, class Compare>
	inline PairingHeap<Key, Compare>::PairingHeap(const Compare &function) : root(0), nodes(0), function(function) { }

	template <class Key, class Compare>
	PairingHeap<Key, Compare>::~PairingHeap()
	{
		clear();
	}

	template <class Key, class Compare>
	typename PairingHeap<Key, Compare>::value_type& PairingHeap<Key, Compare>::top() const
	{
		return root->key;
	}

	template <class Key, class Compare>
	typename PairingHeap<Key, Compare>::iterator PairingHeap<Key, Compare>::push(const typename PairingHeap<Key, Compare>::value_type &key)
	{
		nodes++;
		node *A = new node(key);

		if(root == 0)
			return iterator(root = A);

		if(!function(A->key, root->key))
			A->insert(root), root = A;
		else
			root->insert(A);
		return iterator(A);
	}

	template <class Key, class Compare>
	void PairingHeap<Key, Compare>::pop()
	{
		nodes--;

		if(nodes == 0)
		{
			delete root, root = 0;
			return;
		}

		node *A = root->child, *B, *C;
		delete root, root = A, root->parent = 0;

		while(A)
		{
			B = A->next;
			if(!B)
				break;

			C = B->next;
			if(!function(A->key, B->key))
			{
				if(B->next)
					B->next->previous = A;
				A->next = B->next, A->insert(B), A = A->next;
			}
			else
			{
				if(A->previous)
					A->previous->next = B;
				B->previous = A->previous, B->insert(A), A = B->next;
			}
		}

		if(root->parent)
			root = root->parent;
		A = root->next;
		while(A)
		{
			if(!function(A->key, root->key))
				A->insert(root), A->previous = 0, root = A;
			else
				root->next = A->next, root->insert(A);

			A = root->next;
		}
		root->parent = 0;
	}

	template <class Key, class Compare>
	void PairingHeap<Key, Compare>::update(const typename PairingHeap<Key, Compare>::iterator element, const typename PairingHeap<Key, Compare>::value_type &key)
	{
		node *A = element.data;
		assert(!function(key, A->key));

		A->key = key;

		if(!A->parent)
			return;
		A->remove();

		if(!function(A->key, root->key))
			A->insert(root), root = A;
		else
			root->insert(A);
	}

	template <class Key, class Compare>
	void PairingHeap<Key, Compare>::erase(const typename PairingHeap<Key, Compare>::iterator element)
	{
		node *A = element.data;
		if(A->parent)
			A->remove(), A->insert(root), A->parent = 0, root = A;
		pop();
	}

	template <class Key, class Compare>
	void PairingHeap<Key, Compare>::clear()
	{
		if(root)
			root->clear(), root = 0, nodes = 0;
	}

	template <class Key, class Compare>
	unsigned PairingHeap<Key, Compare>::size() const
	{
		return nodes;
	}

	template <class Key, class Compare>
	bool PairingHeap<Key, Compare>::empty() const
	{
		return !root;
	}
}

#endif
