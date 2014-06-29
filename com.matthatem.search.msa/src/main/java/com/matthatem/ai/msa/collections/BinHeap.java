/**
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package com.matthatem.ai.msa.collections;

import java.util.ArrayList;
import java.util.Comparator;

/**
 * An implementation of a binary heap where elements are aware of their 
 * location (index) in the heap.
 * 
 * @author Matthew Hatem
 */
public class BinHeap<T extends Indexable> {
	
  private final ArrayList<T> heap;
	private final Comparator<T> cmp;
	
	public BinHeap(Comparator<T> cmp) {
	  this.heap = new ArrayList<T>();
		this.cmp = cmp;
	}
	
	public boolean isEmpty() {
		return heap.isEmpty();
	}	
	
	public int size() {
	  return heap.size();
	}
	
	public T poll() {
		if (heap.isEmpty())
		  return null;		
		T t = heap.get(0);
    setIndex(t, -1);    
		if (heap.size() > 1) {
		  T b = heap.remove(heap.size()-1);
		  heap.set(0, b);
		  setIndex(b, 0);
		  pushDown(0);
		}
		else {
		  heap.remove(0);
		}
		return t;
	}
	
	public void add(T i) {
	  heap.add(i);
	  setIndex(i, heap.size()-1);
	  pullUp(heap.size()-1);
	}
	
	public void clear() {
	  heap.clear();
	}
	
	public void update(int i) {
	  if (i < 0 || i > heap.size())
	    throw new IllegalArgumentException();
	  i = pullUp(i);
	  pushDown(i);
	}	
	
	private int pullUp(int i) {
		if (i == 0) 
		  return i;				
		int p = parent(i);		
		if (compare(i, p) < 0) {
		  swap(i, p);
		  return pullUp(p);
		}
		return i;
	}
	
	private void pushDown(int i) {
	  int l = left(i); 
	  int r = right(i);
	  int sml = i;	  
	  if (l < heap.size() && compare(l, i) < 0)
	    sml = l;
	  if (r < heap.size() && compare(r, sml) < 0)
	    sml = r;
	  if (sml != i) {
	    swap (i, sml);
	    pushDown(sml);
	  }
	}
	
	private int compare(int i, int j) {
	  T a = heap.get(i);
	  T b = heap.get(j);
	  return cmp.compare(a, b);
	}
	
	private void setIndex(T t, int i) {
	  t.setIndex(i);
	}
	
	private void swap(int i, int j) {	  
	  T it = heap.get(i);
	  T jt = heap.get(j);
	  	  
	  heap.set(i, jt);
	  setIndex(jt, i);
	  heap.set(j, it);
	  setIndex(it, j);
	}
	
	private int parent(int i) {
		return (i-1)/2;
	}
	
	private int left (int i) {
	  return 2 * i + 1;
	}
	
	private int right (int i) {
	  return 2 * i + 2;
	}
}
