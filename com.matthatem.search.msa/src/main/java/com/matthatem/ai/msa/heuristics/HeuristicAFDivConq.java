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
package com.matthatem.ai.msa.heuristics;

import com.matthatem.ai.msa.MSAHeuristic;
import com.matthatem.ai.msa.SubMatrix;
import com.matthatem.ai.msa.MSA.MSAState;

/**
 * A heuristic that uses a combination of 2D and 3D heuristics.
 * 
 * @author Matthew Hatem
 */
public class HeuristicAFDivConq implements MSAHeuristic {

  private MSAHeuristic[] h3dTable; 
  private MSAHeuristic[] h2dTable; 
  private int[][] h3dIndex;
  private int[][] h2dIndex;
  
  public static enum TYPE {DOUBLE, INT};
  
  public HeuristicAFDivConq(char[][] seqs, SubMatrix sm, double weight) {
    this(seqs, sm, TYPE.DOUBLE, weight);  
  }
  
  public HeuristicAFDivConq(char[][] seqs, SubMatrix sm, TYPE type, double weight) {
    // build the 3D table
    int n = seqs.length/3;
    if (type == TYPE.INT)
      h3dTable = new HeuristicAF3DInt[n];
    else 
      h3dTable = new HeuristicAF3D[n];
    h3dIndex = new int[n][3];
    
    int h2Size = (3)*(seqs.length-3);
    if (seqs.length>3 && seqs.length%3 == 2) h2Size++;
    h2dTable = new HeuristicAF2D[h2Size]; 
    h2dIndex = new int[h2Size][2];
    
    int h3di=0; int h2di=0;
    for (int i=0; i<n; i++) {
      char[][] seqSub1 = new char[3][];
      h3dIndex[h3di][0] = i*3; 
      h3dIndex[h3di][1] = (i*3)+1; 
      h3dIndex[h3di][2] = (i*3)+2;
      seqSub1[0] = new String(seqs[h3dIndex[h3di][0]]).trim().toCharArray();
      seqSub1[1] = new String(seqs[h3dIndex[h3di][1]]).trim().toCharArray();
      seqSub1[2] = new String(seqs[h3dIndex[h3di][2]]).trim().toCharArray();
      
      if (type == TYPE.INT)
        h3dTable[h3di] = new HeuristicAF3DInt(seqSub1, sm, weight);
      else
        h3dTable[h3di] = new HeuristicAF3D(seqSub1, sm, weight);
      
      h3di++;
      // build the 2D table
      for (int p=0; p<3; p++) {
        for (int j=(i+1)*3; j<seqs.length; j++) {
          char[][] seqSub2 = new char[2][];
          int s = (i*3)+p;
          seqSub2[0] = new String(seqs[s]).trim().toCharArray();
          seqSub2[1] = new String(seqs[j]).trim().toCharArray();
          h2dTable[h2di] = new HeuristicAF2D(seqSub2, sm, weight);
          h2dIndex[h2di][0] = s; h2dIndex[h2di][1] = j;
          h2di++;
        }
      }
    }
    
    // if the last set has two sequences
    if (seqs.length>3 && seqs.length%3 == 2) {
      char[][] seqSub2 = new char[2][];
      int s = seqs.length-1; int j = seqs.length-2;
      seqSub2[0] = new String(seqs[s]).trim().toCharArray();
      seqSub2[1] = new String(seqs[j]).trim().toCharArray();
      h2dTable[h2di] = new HeuristicAF2D(seqSub2, sm, weight);
      h2dIndex[h2di][0] = s; h2dIndex[h2di][1] = j;
    }
  }
  
  public double getInitH() {
    double initH = 0;
    for (int i=0; i<h3dTable.length; i++)
      initH += h3dTable[i].getInitH();
    for (int i=0; i<h2dTable.length; i++)
      initH += h2dTable[i].getInitH();
    return initH;
  }

  public double getH(MSAState state, int[] delta, int[] index) {
    throw new IllegalArgumentException();
  }
  
  public double getH(MSAState state, int[] delta) {
    double h = 0;
    for (int i=0; i<h3dTable.length; i++)
      h += h3dTable[i].getH(state, delta, h3dIndex[i]);
    for (int i=0; i<h2dTable.length; i++)
      h += h2dTable[i].getH(state, delta, h2dIndex[i]);
    return h;
  }

}
