#pragma once

//C++简易大顶堆
//Author:Qichao Ying


#ifndef maxheap_h
#define maxheap_h
template<class T>
class Maxheap
{
public:
	Maxheap(int size);
	~Maxheap();
	bool Isempty();
	void push(T item);  //插入操作
	void pop();  //删除操作
	void clear();  //删除操作
	T top();
	T *heap;
	int currentSize;
	int capacity;
};
//-------构造函数初始化-------
template<class T>
Maxheap<T>::Maxheap(int size)
{
	if (size < 1)
	{
		throw"capacity must be >=1";
	}
	else
	{
		currentSize = 0;
		capacity = size;
		heap = new T[capacity + 1]; //heap[0]不使用
	}
}
//-------析构函数释放空间-------
template<class T>
Maxheap<T>::~Maxheap()
{
	delete[]heap;
}
//--------判断堆是否为空-------
template<class T>
bool Maxheap<T>::Isempty()
{
	return currentSize == 0;
}
//---------获取最大元素----------
template<class T>
T Maxheap<T>::top()
{
	return heap[1];
}
//-------插入操作-----
template<class T>
void Maxheap<T>::push(T item)
{
	if (currentSize == capacity)
		throw "Maxheap is full";
	else
	{
		currentSize++;
		int currentNode = currentSize;// 元素的插入位置初始化为最后
		while (currentNode > 1 && heap[currentNode / 2] < item)  //(从下往上)进行调整
		{
			heap[currentNode] = heap[currentNode / 2];
			currentNode = currentNode / 2;
		}
		heap[currentNode] = item; //插入元素
	}
}

template<class T>
void Maxheap<T>::clear()
{
	heap = new T[capacity + 1];
	currentSize=0;
}

//-----删除操作-------
template<class T>
void Maxheap<T>::pop()
{
	if (Isempty())
		throw"heap is empty ,cannot delete";
	else
	{
		T last = heap[currentSize];  //将最后一个元素初始化为根
		currentSize--;
		int currentNode = 1;
		int child = 2;
		while (child < currentSize)  //（从上往下）进行调整
		{
			if (child < currentSize&&heap[child] < heap[child + 1])
				child++;
			if (last > heap[child])
				break;
			else
			{
				heap[currentNode] = heap[child];
				currentNode = child;
				child = child * 2;
			}
		}
		heap[currentNode] = last;
	}
}
#endif

//Example:
//Maxheap<entity> H1(3);
//for (int i = 0; i < 3; i++) {
//	entity e(i, -i);
//	H1.push(e);
//}
//
//cout << endl << H1.top().value << endl;

