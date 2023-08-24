#include <cassert>
#include <cstring>
#include <functional>

namespace Koala
{
	template <class Key, class Compare = std::less<Key> >
	class FibonacciHeap
	{
	public:
		class node
		{
			friend class FibonacciHeap<Key, Compare>;
		
			node *parent, *child, *previous, *next;
			unsigned flag;
		public:
			Key key;
		private:	
			inline node(const Key& = Key());
			
			inline void insert(node*);
			inline void remove();
			
			void clear();
			void check() const;
		};
	private:
		node *root, **degrees;
		unsigned nodes;
		Compare function;
	public:
		inline FibonacciHeap(const Compare& = Compare());
		~FibonacciHeap();
			
		Key& top() const;
		node* push(const Key&);
		void pop();
		
		void update(node*, const Key&);
		void erase(node*);
		
		void merge(FibonacciHeap&);
		void clear();
		
		unsigned size() const;
		bool empty() const;
		
		void check() const;
	};

template <class Key, class Compare>
inline FibonacciHeap<Key, Compare>::node::node(const Key &key) : parent(0), child(0), previous(this), next(this), flag(0), key(key) { }

template <class Key, class Compare>
inline void FibonacciHeap<Key, Compare>::node::insert(typename FibonacciHeap<Key, Compare>::node *A)
{
	next->previous = A->previous, A->previous->next = next, next = A, A->previous = this;
}

template <class Key, class Compare>
inline void FibonacciHeap<Key, Compare>::node::remove()
{
	previous->next = next, next->previous = previous, previous = next = this;
}

template <class Key, class Compare>
void FibonacciHeap<Key, Compare>::node::clear()
{
	if(child)
		child->clear();
	if(previous != this)
		next->previous = previous, next->clear();
	
	delete this;
}

template <class Key, class Compare>
void FibonacciHeap<Key, Compare>::node::check() const
{
	node *child = this->child;
	unsigned degree = 0;

	if(child)
		do
		{
			if(child->previous->next != child || child->next->previous != child)
				throw "Invalid child linked list";
			if(child->parent != this)
				throw "Invalid child parent";
			child->check(), child = child->next, degree++;
		}
		while(child != this->child);
	if(degree != (flag >> 1))
		throw "Invalid child information";
}

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
Key& FibonacciHeap<Key, Compare>::top() const
{
	return root->key;
}

template <class Key, class Compare>
typename FibonacciHeap<Key, Compare>::node* FibonacciHeap<Key, Compare>::push(const Key &key)
{
	nodes++;
	
	node *A = new node(key);
	if(!root)
		return root = A;

	root->insert(A);
	if(function(A->key, root->key))
		root = A;
	return A;
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
	
	if(!nodes)
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
			if(function(C->key, A->key))
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
			if(function(degrees[degree]->key, root->key))
				root = degrees[degree];
			degrees[degree] = 0;
		}
}

template <class Key, class Compare>
void FibonacciHeap<Key, Compare>::update(node *A, const Key &key)
{
	assert(function(key, A->key));

	A->key = key;
	node *B = A->parent;
	if(!B)
	{
		if(function(key, root->key))
			root = A;
		return;
	}
	else if(!function(key, B->key))
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
		if(function(A->key, root->key))
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
void FibonacciHeap<Key, Compare>::erase(node *A)
{
	node *B = A->parent, *C = A;
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
void FibonacciHeap<Key, Compare>::merge(FibonacciHeap& heap)
{
	if(!heap.root || root == heap.root)
		return;
	else if(root)
	{
		root->insert(heap.root);
		if(function(heap.root->key, root->key))
			root = heap.root;
		nodes += heap.nodes;
	}
	else
		root = heap.root, nodes = heap.nodes;
	heap.root = 0, heap.nodes = 0;
}

template <class Key, class Compare>
void FibonacciHeap<Key, Compare>::clear()
{
	if(root)
		root->clear();
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

template <class Key, class Compare>
void FibonacciHeap<Key, Compare>::check() const
{
	if(!root)
		return;

	node *A = root;
	do
	{
		if(A->next->previous != A)
			throw "Invalid root linked list";
		if(A->parent)
			throw "Invalid root parent";
		A->check(), A = A->next;
	}
	while(A != root);
}

}
