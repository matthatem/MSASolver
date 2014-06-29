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
package com.matthatem.ai.msa;

import com.matthatem.ai.msa.MSA.MSAState;

/**
 * The MSA lattice projection class.  This class can be used to form an
 * abstraction of the MSA state space.
 * 
 * @author Matthew Hatem
 */
class MSALatticeProjection implements MSAProjection {
  private int abstraction[][];
  private int abstractionSize;
  private int longestSeqIndex;
  private int longestSeqIndex2;
  
  public MSALatticeProjection (int abs, int absMod, int longestSeqLength, int longestSeqLength2, 
      int longestSeqIndex, int longestSeqIndex2) {
    this.longestSeqIndex = longestSeqIndex;
    this.longestSeqIndex2 = longestSeqIndex2;
    if (abs == 2) {
      int x = longestSeqLength+1;
      int y = longestSeqLength2+1;
      abstractionSize = (absMod > 0) ? absMod : (x*y);
      abstraction = new int[x][y];
      int id = 0;
      for (int i = 0; i < x; i++) {
        for (int j = 0; j < y; j++) {
          abstraction[i][j] = (absMod > 0) ? id % absMod : id;
          id++;
        }
      }
    }
    else {
      int x = longestSeqLength+1;
      int y = longestSeqLength2+1;
      abstractionSize = x;
      abstraction = new int[x][y];
      for (int i = 0; i < x; i++) {
        for (int j = 0; j < y; j++) {
          abstraction[i][j] = i;
        }
      }
    }
  }
  
  public long project(MSAState state) {
    return abstraction[state.pos[longestSeqIndex]]
        [state.pos[longestSeqIndex2]];
  }

  public long size() {
    return abstractionSize;
  }

}
