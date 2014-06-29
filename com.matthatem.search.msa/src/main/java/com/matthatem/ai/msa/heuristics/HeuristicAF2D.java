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
 * A 2D heuristic for MSA with affine gap costs.
 * 
 * @author Matthew Hatem
 */
public class HeuristicAF2D implements MSAHeuristic {
   
  private static final int VT = 0;
  private static final int HZ = 1;
  private static final int DG = 2;
  
  private double afGapCost;
  private double tmGapCost;
  private double lrGapCost;
  private double weight;
  
  private char[][] seqs;
  private double D[][];
  private double H2[][][];
  private double[][][][][] scoreTable;
    
  public HeuristicAF2D(char[][] seqs, SubMatrix sm, double weight) {
    D = sm.getTable();
    this.seqs = seqs;
    lrGapCost = sm.getLinearGapCost();
    afGapCost = sm.getAffineGapCost();
    tmGapCost = sm.getTerminalGapCost();
    this.weight = weight;
    compute();
  }

  private void compute() {
    scoreTable = new double[seqs.length][seqs.length][][][];
    for (int i=1; i<seqs.length; i++) {
      for (int j=0; j<i; j++) {
        char[] A = new String(seqs[i]).trim().toCharArray();
        char[] B = new String(seqs[j]).trim().toCharArray();
        scoreTable[i][j] = DP2(A, B);
      }
    }
    initGapTable();
  }
  
  private void initGapTable() {
    H2 = new double[3][3][3];
    int i, j, k;
    for (i = 0; i < 3; ++i) {
      for (j = 0; j < 3; ++j) {
        for (k = 0; k < 3; ++k) {
          H2[i][j][k] = 0;          
        }
      }
    }       
    H2[0][1][HZ] = afGapCost;   /* two consecutive horizontal moves */
    H2[1][0][VT] = afGapCost;   /* two consecutive vertical moves */
  }
  
  private double[][][] DP2(char[] A, char[] B) {
    int i, j;
    int pos, end, pos2;
    int n = A.length;
    int m = B.length;
    double gapH, gapV;
    double P[][][] = new double [n+1][m+1][3];

    end = m;
    P[n][end][DG] = 0;
    P[n][end][HZ] = P[n][end][VT] = tmGapCost;
    pos = end-1;
    for (j = m - 1; j >= 0; --j, pos -= 1) {
      P[n][pos][VT] = (P[n][pos][DG] = 
                       P[n][pos][HZ] = 
                       P[n][pos+1][HZ] + lrGapCost) + afGapCost;
    }
    for (i = n - 1; i >= 0; --i) {
      gapH = (i == 0 ? tmGapCost : afGapCost);

      P[i][end][HZ] = (P[i][end][DG] = 
                       P[i][end][VT] = 
                       P[i + 1][end][VT] + lrGapCost) + afGapCost;

      pos2 = end;
      pos = end-1;
      for (j = m - 1; j >= 0; --j, pos -= 1) {
        gapV = (j == 0 ? tmGapCost : afGapCost);
        P[i][pos][DG] = min(P[i + 1][pos2][DG], 
                            P[i + 1][pos2][HZ],
                            P[i + 1][pos2][VT]) + D[A[i]][B[j]];

        P[i][pos][HZ] = min(P[i][pos2][DG] + gapH, 
                            P[i][pos2][HZ],
                            P[i][pos2][VT] + gapH) + lrGapCost;

        P[i][pos][VT] = min(P[i + 1][pos][DG] + gapV, 
                            P[i + 1][pos][HZ] + gapV, 
                            P[i + 1][pos][VT]) + lrGapCost;
        pos2 = pos;
      }
    } 
    return P;
  }
  
  public double getH(MSAState state, int[] delta, int[] index) {
    int[] pos = state.pos;
    double cost = 0.0f;
    for (int i=1; i<seqs.length; i++) {
      for (int j = 0; j < i; j++) {
        int di = index[i]; int dj = index[j];
        int col = pos[index[i]]; int row = pos[index[j]];
        cost += min(scoreTable[i][j][col][row][HZ]-H2[delta[di]][delta[dj]][HZ],
                    scoreTable[i][j][col][row][VT]-H2[delta[di]][delta[dj]][VT],
                    scoreTable[i][j][col][row][DG]);
      }
    }
    return weight*cost;
  }
  
  public double getH(MSAState state, int[] delta) {
    int[] pos = state.pos;
    double cost = 0.0f;
    for (int i=1; i<seqs.length; i++) {
      for (int j = 0; j < i; j++) {
        int col = pos[i]; int row = pos[j];
        cost += min(scoreTable[i][j][col][row][HZ]-H2[delta[i]][delta[j]][HZ],
                    scoreTable[i][j][col][row][VT]-H2[delta[i]][delta[j]][VT],
                    scoreTable[i][j][col][row][DG]);
      }
    }
    return weight*cost;
  }
  
  public double getInitH() {
    double cost = 0.0f;
    for (int i=1; i<seqs.length; i++) {
      for (int j = 0; j < i; j++) {
        cost += min(scoreTable[i][j][0][0][HZ],
                    scoreTable[i][j][0][0][VT],
                    scoreTable[i][j][0][0][DG]);
      }
    }
    return weight*cost;
  }
  
  private static final double min(double x, double y, double z) {
    return Math.min(Math.min(x, y), z);
  }
  
  public String toString() {
    StringBuffer sb = new StringBuffer();    
    for (int i=1; i<seqs.length; i++) {
      for (int j=0; j<i; j++) {
        double t[][][] = scoreTable[i][j];
        for (int x=0; x<t.length; x++) {
          for (int y=0; y<t[x].length; y++) {
            sb.append("["+(int)t[x][y][DG]+"|");
            sb.append((int)t[x][y][HZ]+"|");
            sb.append((int)t[x][y][VT]+"]  ");
          }
          sb.append("\n");
        }
        sb.append("\n\n\n");
      }
    }    
    return sb.toString();
  }

}
