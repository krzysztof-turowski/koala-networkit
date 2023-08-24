/*
 * BinomialHeap.hpp
 *
 *  Created on: 24.08.2023
 *      Author: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#pragma once

#include <cassert>
#include <functional>

namespace Koala
{
	template <class Key, class Compare = std::less<Key> >
	class BinomialHeap
	{
	public:
		class node
		{
			friend class BinomialHeap<Key, Compare>;
		
			node *parent, *child, *next;
			unsigned degree;
		public:
			Key key;
		private:	
			inline node(const Key& = Key());
			
			inline void insert(node*);
			
			void clear();
			void check() const;
		};
	private:
		node *root, *minimum;
		unsigned nodes;
		Compare function;
	public:
		inline BinomialHeap(const Compare& = Compare());
		~BinomialHeap();
			
		Key& top() const;
		node* push(const Key&);
		void pop();
		
		void update(node*, const Key&);
		void erase(node*);
		
		void merge(BinomialHeap&);
		void clear();
		
		unsigned size() const;
		bool empty() const;
		
		void check() const;
	private:
		inline node* join(node*, node*);
		inline node* reverse(node*);
		inline node* cut(node*);
	};

template <class Key, class Compare>
inline BinomialHeap<Key, Compare>::node::node(const Key &key) : parent(0), child(0), next(0), degree(0), key(key) { }

template <class Key, class Compare>
inline void BinomialHeap<Key, Compare>::node::insert(typename BinomialHeap<Key, Compare>::node *A)
{
	A->parent = this, A->next = child, child = A, degree++;
}

template <class Key, class Compare>
void BinomialHeap<Key, Compare>::node::clear()
{
	if(next)
		next->clear();
	if(child)
		child->clear();
	
	delete this;
}

template <class Key, class Compare>
void BinomialHeap<Key, Compare>::node::check() const
{
	unsigned degree = 0;
	for(node *child = this->child; child; child = child->next, degree++)
	{
		if(child->parent != this)
			throw "Invalid child parent";
		child->check();
	}
	if(degree != this->degree)
		throw "Invalid child information";
}

template <class Key, class Compare>
inline BinomialHeap<Key, Compare>::BinomialHeap(const Compare &function) : root(0), minimum(0), nodes(0), function(function) { }

template <class Key, class Compare>
BinomialHeap<Key, Compare>::~BinomialHeap()
{
	clear();
}

template <class Key, class Compare>
Key& BinomialHeap<Key, Compare>::top() const
{
	return minimum->key;
}

template <class Key, class Compare>
typename BinomialHeap<Key, Compare>::node* BinomialHeap<Key, Compare>::push(const Key &key)
{
	nodes++;
	node *A = new node(key);

	if(root == 0)
		return root = minimum = A;

	root = join(root, A);
	if(function(A->key, minimum->key))
		minimum = A;
	return A;
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
			if(function(A->key, minimum->key))
				minimum = A;
}

template <class Key, class Compare>
void BinomialHeap<Key, Compare>::update(node *A, const Key &key)
{
	assert(function(key, A->key));
	
	A->key = key;
	if(function(key, minimum->key))
		minimum = A;

	if(!A->parent || function(A->parent->key, A->key))
		return;

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
}

template <class Key, class Compare>
void BinomialHeap<Key, Compare>::erase(node *A)
{
	nodes--;
	
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
			if(function(B->key, minimum->key))
				minimum = B;
	}
	delete A;
}

template <class Key, class Compare>
void BinomialHeap<Key, Compare>::merge(BinomialHeap& heap)
{
	if(!heap.root || root == heap.root)
		return;
	else if(root)
	{
		root = join(root, heap.root);
		if(function(heap.minimum->key, minimum->key))
			minimum = heap.minimum;
		nodes += heap.nodes;
	}
	else
		root = heap.root, minimum = heap.minimum, nodes = heap.nodes;
	heap.root = heap.minimum = 0, heap.nodes = 0;
}

template <class Key, class Compare>
void BinomialHeap<Key, Compare>::clear()
{
	if(root)
		root->clear();
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
void BinomialHeap<Key, Compare>::check() const
{
	node *A = root;
	while(A)
	{
		if(A->parent)
			throw "Invalid root parent";
		A->check(), A = A->next;
	}
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
		else if(function(B->key, C->key))
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
