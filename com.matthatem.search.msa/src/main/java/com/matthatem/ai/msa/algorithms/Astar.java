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
package com.matthatem.ai.msa.algorithms;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

import com.carrotsearch.hppc.LongObjectOpenHashMap;
import com.matthatem.ai.msa.MSA;
import com.matthatem.ai.msa.MSA.MSANode;
import com.matthatem.ai.msa.MSA.MSAState;
import com.matthatem.ai.msa.collections.BinHeap;
import com.matthatem.ai.msa.collections.Indexable;


/**
 * An implementation of A* tailored to the MSA domain.
 * 
 * @author Matthew Hatem
 */
public final class Astar implements SearchAlgorithm {
  
  private LongObjectOpenHashMap closed[];
  private BinHeap<MNH> open = new BinHeap<MNH>(new NodeComparator());
  private List<MSAState> path = new ArrayList<MSAState>(3);
  private MSA domain;
  private long expanded;
  private long generated;
  private long duplicates;
  private MSAState state;
  private MSANode children[];
  
  /**
   * The constructor.
   * 
   * @param domain the search domain
   */
  public Astar(final MSA domain) {
    this.domain = domain;
    this.children = new MSANode[domain.getNumOps()];
    this.closed = new LongObjectOpenHashMap[domain.getNumOps()+1];
    this.closed[domain.getNumOps()] = new LongObjectOpenHashMap<MNH>(); // root op
    for (int i=0; i<this.children.length; i++) {
      this.children[i] = new MSANode();
      this.closed[i] = new LongObjectOpenHashMap<MNH>();
    }
  }
  
  /* (non-Javadoc)
   * @see edu.unh.ai.search.SearchAlgorithm#search(java.lang.Object)
   */
  public SearchResult<MSAState> search() {
    double cost = 0;
    state = domain.initial();    
    MSANode initNode = domain.initialNode();    
    open(new MNH(initNode));
    
    while (!open.isEmpty() && path.isEmpty()) {  
      MNH mnh = open.poll();
      MSANode n = mnh.node;
      if (n == domain.goal) {
        cost = n.g;
        for (MSANode p = n; p != null; p = p.parent) {
          path.add(domain.unpack(p.packed, p.e));
        }
        break;
      }
      expanded++;
      
      domain.unpack(n.packed, n.e, state);
      int count = domain.expand(state, n, children);
      
      for (int i=0; i<count; i++) {        
        MSANode child = children[i];
        children[i] = new MSANode(); // reset!
        // merge duplicates
        MNH dup = (MNH)closed[child.e].get(child.packed);
        if (closed[child.e].containsKey(child.packed)) {
          duplicates++;
          if (child.g >= dup.node.g) {
            continue;
          }
          else {
            dup.node.f = child.f;
            dup.node.g = child.g;
            dup.node.parent = child.parent;
            if (dup.index != -1) {
              open.update(dup.index);
            }
            else {
              open.add(dup);
            }
          }
        }
        // no duplicates
        else {
          open(new MNH(child));
          generated++;
        }        
      }
    }

    SearchResult<MSAState> result = 
        new SearchResult<MSAState>(path, expanded, generated);
    result.setProperty("duplicates", Long.toString(duplicates));
    result.setCost(cost);
    return result;
  }
  
  private void open(MNH mnh) {
    open.add(mnh);
    closed[mnh.node.e].put(mnh.node.packed, mnh);
  }
  
  /*
   * Wrapper for MSANode to make heapable
   */
  private final class MNH implements Indexable {
    private MSANode node;
    private int index = -1;
    public MNH(MSANode node) {
      this.node = node;
    }
    public int getIndex() {
      return index;
    }
    public void setIndex(int index) {
      this.index = index;
    }    
  }
  
  /*
   * The node comparator class
   */
  private final class NodeComparator implements Comparator<MNH> {
    public int compare(final MNH a, final MNH b) {
      if (a.node.f == b.node.f) { 
        if (a.node.g > b.node.g) return -1;
        if (a.node.g < b.node.g) return 1;
        return 0;
      }
      else {
        if (a.node.f < b.node.f) return -1;
        if (a.node.f > b.node.f) return 1;
        return 0;
      }
    }    
  }
  
}
