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
 * A 3D heuristic for MSA with affine gap costs and integer values.
 * 
 * @author Matthew Hatem
 */
public class HeuristicAF3DInt implements MSAHeuristic {
    
  private static final int X = 0;
  private static final int Y = 1;
  private static final int Z = 2;
  private static final int XY = 3;
  private static final int YZ = 4;
  private static final int ZX = 5;
  private static final int XYZ = 6;
  
  private int afGapCost;
  private int tmGapCost;
  private int lrGapCost;
  private double weight;
  
  private char[][] seqs;
  private double D[][];
  private int H2[][][][];
  private int[][][][][][][] scoreTable;
  
  private int n;
  private int[][] ni = new int[3][3];
  
  
  public HeuristicAF3DInt(char[][] seqs, SubMatrix sm, double weight) {
    D = sm.getTable();
    this.seqs = seqs;
    lrGapCost = (int)sm.getLinearGapCost();
    afGapCost = tmGapCost = (int)sm.getAffineGapCost();
    this.weight = weight;
    compute();
  }

  private void compute() {
    scoreTable = new int[seqs.length][seqs.length][seqs.length][][][][];
    n = seqs.length/3;
    for (int i=0; i<n; i++) {
      ni[i][0] = i*3; ni[i][1] = (i*3)+1; ni[i][2] = (i*3)+2;
      char[] A = new String(seqs[ni[i][0]]).trim().toCharArray();
      char[] B = new String(seqs[ni[i][1]]).trim().toCharArray();
      char[] C = new String(seqs[ni[i][2]]).trim().toCharArray();
      scoreTable[ni[i][0]][ni[i][1]][ni[i][2]] = DP2(A, B, C);
    }
    initGapTable();
  }
  
  private void initGapTable() {
    H2 = new int[2][2][2][7];
    int h, i, j, k;
    for (h = 0; h < 2; ++h) {
      for (i = 0; i < 2; ++i) {
        for (j = 0; j < 2; ++j) {
          for (k = 0; k < 7; ++k) {
            H2[h][i][j][k] = 0;
          }
        }
      }
    }
    
    //[X][Y][Z]
    H2[1][0][0][X] = afGapCost*2;
    H2[1][1][0][X] = afGapCost; 
    H2[1][0][1][X] = afGapCost; 
    
    H2[0][1][0][Y] = afGapCost*2;   
    H2[1][1][0][Y] = afGapCost;
    H2[0][1][1][Y] = afGapCost;
    
    H2[0][0][1][Z] = afGapCost*2;   
    H2[1][0][1][Z] = afGapCost;
    H2[0][1][1][Z] = afGapCost;
    
    H2[1][0][0][XY] = afGapCost;   
    H2[0][1][0][XY] = afGapCost;
    H2[1][1][0][XY] = afGapCost*2;
    
    H2[0][1][0][YZ] = afGapCost;   
    H2[0][0][1][YZ] = afGapCost;
    H2[0][1][1][YZ] = afGapCost*2;
    
    H2[1][0][0][ZX] = afGapCost;
    H2[0][0][1][ZX] = afGapCost;
    H2[1][0][1][ZX] = afGapCost*2;
  }
  
  private int[][][][] DP2(final char[] A, final char[] B, final char[] C) {    
    final int x = A.length;
    final int y = B.length;
    final int z = C.length;
    int P[][][][] = new int [x+1][y+1][z+1][7];
    
    // initialize the starting cell
    P[x][y][z][X] = P[x][y][z][Y] = P[x][y][z][Z] = 2*tmGapCost;
    P[x][y][z][XY] = P[x][y][z][YZ] = P[x][y][z][ZX] = 2*tmGapCost;
    P[x][y][z][XYZ] = 0;
    
    // initialize the x edge
    for (int xpos = x-1; xpos >= 0; --xpos) {
      P[xpos][y][z][X] = P[xpos+1][y][z][X]+2*lrGapCost;
      P[xpos][y][z][Y] = P[xpos+1][y][z][X]+2*lrGapCost+2*afGapCost;
      P[xpos][y][z][Z] = P[xpos+1][y][z][X]+2*lrGapCost+2*afGapCost;
      P[xpos][y][z][XY] = P[xpos+1][y][z][X]+2*lrGapCost+afGapCost;
      P[xpos][y][z][YZ] = P[xpos+1][y][z][X]+2*lrGapCost+2*afGapCost;
      P[xpos][y][z][ZX] = P[xpos+1][y][z][X]+2*lrGapCost+afGapCost;
      P[xpos][y][z][XYZ] = P[xpos+1][y][z][X]+2*lrGapCost;
      
    }
       
    // initialize the y edge
    for (int ypos = y-1; ypos >= 0; --ypos) {
      P[x][ypos][z][X] = P[x][ypos+1][z][Y]+2*lrGapCost+2*afGapCost;
      P[x][ypos][z][Y] = P[x][ypos+1][z][Y]+2*lrGapCost;
      P[x][ypos][z][Z] = P[x][ypos+1][z][Y]+2*lrGapCost+2*afGapCost;
      P[x][ypos][z][XY] = P[x][ypos+1][z][Y]+2*lrGapCost+afGapCost;
      P[x][ypos][z][YZ] = P[x][ypos+1][z][Y]+2*lrGapCost+afGapCost;
      P[x][ypos][z][ZX] = P[x][ypos+1][z][Y]+2*lrGapCost+2*afGapCost;
      P[x][ypos][z][XYZ] = P[x][ypos+1][z][Y]+2*lrGapCost;
    }
    
    // initialize the z edge
    for (int zpos = z-1; zpos >= 0; --zpos) {
      P[x][y][zpos][X] = P[x][y][zpos+1][Z]+2*lrGapCost+2*afGapCost;
      P[x][y][zpos][Y] = P[x][y][zpos+1][Z]+2*lrGapCost+2*afGapCost;
      P[x][y][zpos][Z] = P[x][y][zpos+1][Z]+2*lrGapCost;
      P[x][y][zpos][XY] = P[x][y][zpos+1][Z]+2*lrGapCost+2*afGapCost;
      P[x][y][zpos][YZ] = P[x][y][zpos+1][Z]+2*lrGapCost+afGapCost;
      P[x][y][zpos][ZX] = P[x][y][zpos+1][Z]+2*lrGapCost+afGapCost;
      P[x][y][zpos][XYZ] = P[x][y][zpos+1][Z]+2*lrGapCost;
    }
    
    // init the xy face
    for (int xpos = x-1; xpos >= 0; --xpos) {
      for (int ypos = y-1; ypos >= 0; --ypos) {
        // x
        P[xpos][ypos][z][X] = 
           min(P[xpos+1][ypos][z][XY]+(afGapCost),
               P[xpos+1][ypos][z][X],
               P[xpos+1][ypos][z][Y]+(afGapCost))
           +(2*lrGapCost);        
        // y
        P[xpos][ypos][z][Y] = 
            min(P[xpos][ypos+1][z][XY]+(afGapCost),
                P[xpos][ypos+1][z][X]+(afGapCost),
                P[xpos][ypos+1][z][Y])
            +(2*lrGapCost);
        // xy
        P[xpos][ypos][z][XY] = 
            min(P[xpos+1][ypos+1][z][XY],
                P[xpos+1][ypos+1][z][X]+(afGapCost),
                P[xpos+1][ypos+1][z][Y]+(afGapCost))
            +(2*lrGapCost)+
            (int)D[A[xpos]][B[ypos]];
      
        // this should be Double.MAX_VALUE
        int min = min(P[xpos][ypos][z][X], 
            P[xpos][ypos][z][Y], P[xpos][ypos][z][XY]);
        P[xpos][ypos][z][Z] =  min+2*afGapCost;
        P[xpos][ypos][z][YZ] = min+2*afGapCost;
        P[xpos][ypos][z][ZX] = min+2*afGapCost;
        P[xpos][ypos][z][XYZ]= min+afGapCost;
        
      }
    }
        
    // init the yz face
    for (int ypos = y-1; ypos >= 0; --ypos) {
      for (int zpos = z-1; zpos >= 0; --zpos) {
        // y
        P[x][ypos][zpos][Y] = 
            min(P[x][ypos+1][zpos][YZ]+(afGapCost),
                P[x][ypos+1][zpos][Y],
                P[x][ypos+1][zpos][Z]+(afGapCost))
            +(2*lrGapCost);
        // z
        P[x][ypos][zpos][Z] = 
            min(P[x][ypos][zpos+1][YZ]+(afGapCost),
                P[x][ypos][zpos+1][Y]+(afGapCost),
                P[x][ypos][zpos+1][Z])
            +(2*lrGapCost);
        // yz
        P[x][ypos][zpos][YZ] = 
            min(P[x][ypos+1][zpos+1][YZ],
                P[x][ypos+1][zpos+1][Y]+(afGapCost),
                P[x][ypos+1][zpos+1][Z]+(afGapCost))
            +(2*lrGapCost)+
            (int)D[B[ypos]][C[zpos]];

        // this should be Double.MAX_VALUE
        int min = min(P[x][ypos][zpos][Y], 
            P[x][ypos][zpos][Z], P[x][ypos][zpos][YZ]);
        P[x][ypos][zpos][X] =  min+2*afGapCost;
        P[x][ypos][zpos][XY] = min+2*afGapCost;
        P[x][ypos][zpos][ZX] = min+2*afGapCost;
        P[x][ypos][zpos][XYZ]= min+2*afGapCost;
      }
    }
    
    // init the zx face
    for (int zpos = z-1; zpos >= 0; --zpos) {
      for (int xpos = x-1; xpos >= 0; --xpos) {
        // z
        P[xpos][y][zpos][Z] = 
            min(P[xpos][y][zpos+1][ZX]+(afGapCost),
                P[xpos][y][zpos+1][Z],
                P[xpos][y][zpos+1][X]+(afGapCost))
            +(2*lrGapCost);
        // x
        P[xpos][y][zpos][X] = 
            min(P[xpos+1][y][zpos][ZX]+(afGapCost),
                P[xpos+1][y][zpos][Z]+(afGapCost),
                P[xpos+1][y][zpos][X])
            +(2*lrGapCost);
        // zx
        P[xpos][y][zpos][ZX] = 
            min(P[xpos+1][y][zpos+1][ZX],
                P[xpos+1][y][zpos+1][Z]+(afGapCost),
                P[xpos+1][y][zpos+1][X]+(afGapCost))
            +(2*lrGapCost)+
            (int)D[C[zpos]][A[xpos]];

       // this should be Double.MAX_VALUE
        int min = min(P[xpos][y][zpos][Z], 
            P[xpos][y][zpos][X], P[xpos][y][zpos][ZX]);
        P[xpos][y][zpos][Y] =  min+2*afGapCost;
        P[xpos][y][zpos][XY] = min+2*afGapCost;
        P[xpos][y][zpos][YZ] = min+2*afGapCost;
        P[xpos][y][zpos][XYZ]= min+2*afGapCost;
        
      }
    }
        
    for (int xpos = x - 1; xpos >= 0; --xpos) {
      for (int ypos = y - 1; ypos >= 0; --ypos) {
        for (int zpos = z - 1; zpos >= 0; --zpos) {
          // XYZ
          P[xpos][ypos][zpos][XYZ] = 
              min(P[xpos+1][ypos+1][zpos+1][X], 
                  P[xpos+1][ypos+1][zpos+1][Y],
                  P[xpos+1][ypos+1][zpos+1][Z],
                  P[xpos+1][ypos+1][zpos+1][XY],
                  P[xpos+1][ypos+1][zpos+1][YZ],
                  P[xpos+1][ypos+1][zpos+1][ZX],
                  P[xpos+1][ypos+1][zpos+1][XYZ]) + 
                  (int)D[A[xpos]][B[ypos]]+
                  (int)D[B[ypos]][C[zpos]]+
                  (int)D[C[zpos]][A[xpos]];
          // X
          P[xpos][ypos][zpos][X] = 
              min(P[xpos+1][ypos][zpos][X], 
                  P[xpos+1][ypos][zpos][Y]+afGapCost*2,
                  P[xpos+1][ypos][zpos][Z]+afGapCost*2,
                  P[xpos+1][ypos][zpos][XY]+afGapCost,
                  P[xpos+1][ypos][zpos][YZ]+afGapCost*2,
                  P[xpos+1][ypos][zpos][ZX]+afGapCost,
                  P[xpos+1][ypos][zpos][XYZ]+afGapCost*2) + (2*lrGapCost);
          // Y
          P[xpos][ypos][zpos][Y] = 
              min(P[xpos][ypos+1][zpos][X]+afGapCost*2, 
                  P[xpos][ypos+1][zpos][Y],
                  P[xpos][ypos+1][zpos][Z]+afGapCost*2,
                  P[xpos][ypos+1][zpos][XY]+afGapCost,
                  P[xpos][ypos+1][zpos][YZ]+afGapCost,
                  P[xpos][ypos+1][zpos][ZX]+afGapCost*2,
                  P[xpos][ypos+1][zpos][XYZ]+afGapCost*2) + (2*lrGapCost);
          // Z
          P[xpos][ypos][zpos][Z] = 
              min(P[xpos][ypos][zpos+1][X]+afGapCost*2, 
                  P[xpos][ypos][zpos+1][Y]+afGapCost*2,
                  P[xpos][ypos][zpos+1][Z],
                  P[xpos][ypos][zpos+1][XY]+afGapCost*2,
                  P[xpos][ypos][zpos+1][YZ]+afGapCost,
                  P[xpos][ypos][zpos+1][ZX]+afGapCost,
                  P[xpos][ypos][zpos+1][XYZ]+afGapCost*2) + (2*lrGapCost);
          
          // XY
          P[xpos][ypos][zpos][XY] = 
              min(P[xpos+1][ypos+1][zpos][X]+afGapCost,
                  P[xpos+1][ypos+1][zpos][Y]+afGapCost,
                  P[xpos+1][ypos+1][zpos][Z]+afGapCost*2,
                  P[xpos+1][ypos+1][zpos][XY],
                  P[xpos+1][ypos+1][zpos][YZ]+afGapCost*2,
                  P[xpos+1][ypos+1][zpos][ZX]+afGapCost*2,
                  P[xpos+1][ypos+1][zpos][XYZ]+afGapCost*2) + (2*lrGapCost) +
                  (int)D[A[xpos]][B[ypos]];          
          // YZ
          P[xpos][ypos][zpos][YZ] = 
              min(P[xpos][ypos+1][zpos+1][X]+afGapCost*2, 
                  P[xpos][ypos+1][zpos+1][Y]+afGapCost,
                  P[xpos][ypos+1][zpos+1][Z]+afGapCost,
                  P[xpos][ypos+1][zpos+1][XY]+afGapCost*2,
                  P[xpos][ypos+1][zpos+1][YZ],
                  P[xpos][ypos+1][zpos+1][ZX]+afGapCost*2,
                  P[xpos][ypos+1][zpos+1][XYZ]+afGapCost*2) + (2*lrGapCost) +
                  (int)D[B[ypos]][C[zpos]];
          // ZX
          P[xpos][ypos][zpos][ZX] = 
              min(P[xpos+1][ypos][zpos+1][X]+afGapCost, 
                  P[xpos+1][ypos][zpos+1][Y]+afGapCost*2,
                  P[xpos+1][ypos][zpos+1][Z]+afGapCost,
                  P[xpos+1][ypos][zpos+1][XY]+afGapCost*2,
                  P[xpos+1][ypos][zpos+1][YZ]+afGapCost*2,
                  P[xpos+1][ypos][zpos+1][ZX],
                  P[xpos+1][ypos][zpos+1][XYZ]+afGapCost*2) + (2*lrGapCost) +
                  (int)D[C[zpos]][A[xpos]];
        }
      }
    } 
    return P;
  }
  
  public double getH(MSAState state, int[] delta, int index[]) {
    int[] pos = state.pos;
    int cost = 0;
    for (int p=0; p<n; p++) {
      int h = ni[p][0]; int i = ni[p][1]; int j = ni[p][2];
      int x = pos[index[0]]; int y = pos[index[1]]; int z = pos[index[2]];
      int hi = index[0]; int ii = index[1]; int ji = index[2];
      cost += min(
          scoreTable[h][i][j][x][y][z][X]-H2[delta[hi]][delta[ii]][delta[ji]][X],
          scoreTable[h][i][j][x][y][z][Y]-H2[delta[hi]][delta[ii]][delta[ji]][Y],
          scoreTable[h][i][j][x][y][z][Z]-H2[delta[hi]][delta[ii]][delta[ji]][Z],              
          scoreTable[h][i][j][x][y][z][XY]-H2[delta[hi]][delta[ii]][delta[ji]][XY],
          scoreTable[h][i][j][x][y][z][YZ]-H2[delta[hi]][delta[ii]][delta[ji]][YZ],
          scoreTable[h][i][j][x][y][z][ZX]-H2[delta[hi]][delta[ii]][delta[ji]][ZX],
          scoreTable[h][i][j][x][y][z][XYZ]);
    }
    return weight*cost;
  }
  
  public double getH(MSAState state, int[] delta) {
    int[] pos = state.pos;
    int cost = 0;
    for (int p=0; p<n; p++) {
      int h = ni[p][0]; int i = ni[p][1]; int j = ni[p][2];
      int x = pos[h]; int y = pos[i]; int z = pos[j];
      cost += min(
          scoreTable[h][i][j][x][y][z][X]-H2[delta[h]][delta[i]][delta[j]][X],
          scoreTable[h][i][j][x][y][z][Y]-H2[delta[h]][delta[i]][delta[j]][Y],
          scoreTable[h][i][j][x][y][z][Z]-H2[delta[h]][delta[i]][delta[j]][Z],              
          scoreTable[h][i][j][x][y][z][XY]-H2[delta[h]][delta[i]][delta[j]][XY],
          scoreTable[h][i][j][x][y][z][YZ]-H2[delta[h]][delta[i]][delta[j]][YZ],
          scoreTable[h][i][j][x][y][z][ZX]-H2[delta[h]][delta[i]][delta[j]][ZX],
          scoreTable[h][i][j][x][y][z][XYZ]);
    }
    return weight*cost;
  }
  
  public double getInitH() {
    double cost = 0.0f;
    for (int p=0; p<n; p++) {
      int h = ni[p][0]; int i = ni[p][1]; int j = ni[p][2];
      cost += min(scoreTable[h][i][j][0][0][0][X],
                  scoreTable[h][i][j][0][0][0][Y],
                  scoreTable[h][i][j][0][0][0][Z],                      
                  scoreTable[h][i][j][0][0][0][XY],
                  scoreTable[h][i][j][0][0][0][YZ],
                  scoreTable[h][i][j][0][0][0][ZX],
                  scoreTable[h][i][j][0][0][0][XYZ]);
    }
    return weight*cost;
  }
  
  public int[] getInitRow() {
    int h = ni[0][0]; int i = ni[0][1]; int j = ni[0][2];
    return scoreTable[h][i][j][0][0][0];
  }
  
  private static final int min(int x, int y, int z) {
    return Math.min(Math.min(x, y), z);
  }
  
  private static final int min(int x, int y, int z, int xy, 
      int yz, int zx, int xyz) {
    // TODO maybe something more efficient here!!
    return Math.min(Math.min(Math.min(Math.min(Math.min(Math.min(x, y), z), xy), yz), zx), xyz);
  }
  
  public void printTable(double t[][][][], char[] A, char[] B, char[] C) {
      final int x = A.length;
      final int y = B.length;
      final int z = C.length;
      
      // XY face
      StringBuffer sb = new StringBuffer();
      for (int ypos=y; ypos>=0; ypos--) {
        for (int xpos=x; xpos>=0; xpos--) {
          sb.append("XYZ: "+xpos+" "+ypos+" "+z+"   ");
          sb.append("[");
          for (int q=0; q<6; q++) {
            sb.append((int)t[xpos][ypos][z][q]+"|");
          }
          sb.append((int)t[xpos][ypos][z][6]+"]  ");
        }
        sb.append("\n");
      }
      System.out.println(sb.toString());
      
      // YZ face
      sb = new StringBuffer();
      for (int zpos=z; zpos>=0; zpos--) {
        for (int ypos=y; ypos>=0; ypos--) {
          sb.append("XYZ: "+ypos+" "+zpos+" "+z+"   ");
          sb.append("[");
          for (int q=0; q<6; q++) {
            sb.append((int)t[x][ypos][zpos][q]+"|");
          }
          sb.append((int)t[x][ypos][zpos][6]+"]  ");
        }
        sb.append("\n");
      }
      System.out.println(sb.toString());
      
      // ZX face
      sb = new StringBuffer();
      for (int xpos=x; xpos>=0; xpos--) {
        for (int zpos=z; zpos>=0; zpos--) {
          sb.append("XYZ: "+xpos+" "+zpos+" "+z+"   ");
          sb.append("[");
          for (int q=0; q<6; q++) {
            sb.append((int)t[xpos][y][zpos][q]+"|");
          }
          sb.append((int)t[xpos][y][zpos][6]+"]  ");
        }
        sb.append("\n");
      }
      System.out.println(sb.toString());
  }
  public String toString() {
    StringBuffer sb = new StringBuffer();
    for (int p=0; p<n; p++) {
      int h = ni[p][0]; int i = ni[p][1]; int j = ni[p][2];
      int t[][][][] = scoreTable[h][i][j];
      
      for (int x=0; x<t.length; x++) {
        for (int y=0; y<t[x].length; y++) {
          for (int z=0; z<t[x][y].length; z++) {
            sb.append("XYZ: "+x+" "+y+" "+z+"   ");
            sb.append("[");
            for (int q=0; q<6; q++) {
              sb.append((int)t[x][y][z][q]+"|");
            }
            sb.append((int)t[x][y][z][6]+"]  ");
          }
          sb.append("\n");
        }
        sb.append("\n\n");
      }
      
    }
    sb.append("\n\n\n");
    return sb.toString();
  }

}
