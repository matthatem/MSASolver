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

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;

import com.matthatem.ai.msa.heuristics.HeuristicAF2D;
import com.matthatem.ai.msa.heuristics.HeuristicAF3D;
import com.matthatem.ai.msa.heuristics.HeuristicAF3DInt;
import com.matthatem.ai.msa.heuristics.HeuristicAFDivConq;


/**
 * The MSA domain class.
 * 
 * @author Matthew Hatem
 */
public class MSA {
  
  /*
   * We assume 6 or fewer sequences and the length 
   * of a sequence must be less than 1023.
   */
  private static final int MAX_SEQ_LENGTH = 1023;
  private static final int MAX_SEQ_NUM = 6;
  
  private int[] gcBit, eIDTab;
  private boolean[] gcBitInc;
  private HashMap<String, String> mapIDtoSeq = new HashMap<String, String>();
  private HashMap<Integer, String> mapIndextoID = new HashMap<Integer, String>();
  private MSAState initial;
  private char[][] seqs;
  private int[] seqLen;
  private int numSeqs;
  private int numOps;
  private int longestSeqLength;
  private int longestSeqLength2;
  private int longestSeqIndex;
  private int longestSeqIndex2;
  private MSAHeuristic heuristic;
  private SubMatrix subMatrix;
  private double D[][];
 
  private static final int VT = 0;
  private static final int HZ = 1;
  private double T[][][][];
  private double H2[][][];
  
  private int oprDim[];
  private char oprSeq[][][];
  private char oprSeqPw[][][];
  private int oprSeqPos[][];
  private int oprSeqPwPos[][];
  private double oprGapCost[];
  
  private MSAProjection projection;
  
  private boolean penTermGaps;  
  
  public MSANode goal = null;
  
  public static enum HEURISTICS {H2D, H3D, H3D_INT, HDIVCONQ, HDIVCONQ_INT};
  
  /*
   * The MSA state class.
   */
  public static final class MSAState {  
    public int pos[];
    public double h;
    public byte e;
    public MSAState(int dim) {
      pos = new int[dim];
    }
  }
  
  public static class MSANode {
    public double g;
    public double f;
    public long packed;
    public byte e;
    public MSANode parent;
  }
    
  public MSA(InputStream stream, HEURISTICS h) {
    this(stream, h, true, 1);
  }
  public MSA(InputStream stream, HEURISTICS h, boolean penTermGaps, double weight) {
    this.penTermGaps = penTermGaps;
    try {
      BufferedReader reader = new BufferedReader(
          new InputStreamReader(stream));
      readSequences(reader);
      stream.close();
    } catch (FileNotFoundException e) {
      e.printStackTrace();
    } catch (IOException e) {
      e.printStackTrace();
    }
    numSeqs = mapIDtoSeq.keySet().size();

    this.subMatrix = new SubMatrix("msa/pam250.sub", penTermGaps);
    if (h.equals(HEURISTICS.H3D)) {
      this.heuristic = new HeuristicAF3D(seqs, subMatrix, weight);
    }
    else if (h.equals(HEURISTICS.H3D_INT)) {
      this.heuristic = new HeuristicAF3DInt(seqs, subMatrix, weight);
    }
    else if (h.equals(HEURISTICS.HDIVCONQ)) {
      System.out.println("Using divconq heuristic!");
      this.heuristic = new HeuristicAFDivConq(seqs, subMatrix, weight);
    }
    else if (h.equals(HEURISTICS.HDIVCONQ_INT)) {
      System.out.println("Using divconq_int heuristic!");
      this.heuristic = new HeuristicAFDivConq(seqs, subMatrix,
          HeuristicAFDivConq.TYPE.INT, weight);
    }
    else {
      this.heuristic = new HeuristicAF2D(seqs, subMatrix, weight);
    }
    this.D = subMatrix.D;
    
    initGapTable();
    initGrayCodes();
    initOperators();
  }
  
  private void readSequences(BufferedReader reader) throws IOException {
    String key = null;
    while ((key = reader.readLine()) != null) {
      if (key.trim().length() == 0) continue;
      key = key.substring(2, key.length());
      String seq = " "+reader.readLine();
      if (seq.length()-1 > MAX_SEQ_LENGTH) {
        System.err.println("Sequence exceeds maximum length!");
        System.exit(1);
      }
      mapIDtoSeq.put(key, seq);
      if (mapIDtoSeq.values().size() > MAX_SEQ_NUM) {
        System.err.println("Number of sequences exceeds maximum!");
        System.exit(1);
      }
    }
    seqs = new char[mapIDtoSeq.values().size()][];
    seqLen = new int[mapIDtoSeq.values().size()];
    
    // sort the strings
    ArrayList<Map.Entry<String, String>> list =
        new ArrayList<Map.Entry<String, String>>(mapIDtoSeq.entrySet());
    Collections.sort(list, new Comparator<Map.Entry<String, String>>() {
      public int compare(Entry<String, String> o1, Entry<String, String> o2) {
        return o2.getValue().length() - o1.getValue().length();
      }
    });
    
    for (int id=0; id<list.size(); id++) {
      Map.Entry<String, String> entry = list.get(id);
      String k = entry.getKey();
      String v = entry.getValue();
      if (longestSeqLength == 0) {
      //if ((id+1) == list.size()/2) {
        longestSeqLength = v.length();
        longestSeqIndex = id;
      }
      else if (longestSeqLength2 == 0) {
      //else if ((id+1) == (list.size()/2)+1) {
        longestSeqLength2 = v.length();
        longestSeqIndex2 = id;
      }
      seqs[id] = v.toCharArray();
      seqLen[id] = v.length()-1;
      mapIndextoID.put(id, k);
    }
  }
  
  private void initGrayCodes() {
    int k = numSeqs;
    numOps = (1 << k) - 1;
    int gc[] = new int[numOps+1];
    gcBit = new int[numOps];
    eIDTab = new int[numOps];
    gcBitInc = new boolean[numOps];
    for (int n = 0; n <= numOps; n++)
      gc[n] = n ^ (n >> 1);
    for (int n = 0; n < numOps; n++)
      eIDTab[n] = gc[n + 1];
    int i = 0;
    for (int n = 0; n < numOps; n++) {      
      for (int bit = 0, testBit = 1; bit < numSeqs; bit++, testBit <<= 1) {
        if ((testBit & gc[n]) != (testBit & gc[n+1])) {
          gcBit[i] = bit;
          if ((testBit & gc[n]) > 0)                    
            gcBitInc[i++] = false;
          else                    
            gcBitInc[i++] = true;                       
          break;
        }           
      }
    }
  }
  
  private void initOperators() {
    oprDim = new int[numOps];
    oprSeq = new char[numOps][numSeqs][];
    oprSeqPw = new char[numOps][numSeqs][];
    oprSeqPos = new int[numOps][numSeqs];
    oprSeqPwPos = new int[numOps][numSeqs];
    oprGapCost = new double[numOps];
    
    for (int opr = 0; opr < numOps; ++opr) {
      int n = opr + 1;
      int gc = n ^ (n >> 1);
      oprDim[opr] = 0;
      int pwPosIdx = 0;
      for (int i = 0, testBit = 1; i < numSeqs; ++i, testBit <<= 1) {     
        if ((testBit & gc) > 0) {
          oprSeq[opr][oprDim[opr]] = seqs[i];
          oprSeqPos[opr][oprDim[opr]] = i;
          oprDim[opr]++;
          if (gcBit[opr] != i) {
            oprSeqPw[opr][pwPosIdx] = seqs[i];
            oprSeqPwPos[opr][pwPosIdx] = i;
            pwPosIdx++;
          }
        }
      }
      oprGapCost[opr] = 
        (numSeqs - oprDim[opr]) * oprDim[opr] * subMatrix.getLinearGapCost();
    } 
  }
  
  private void initGapTable() {
    T = new double[3][3][3][3];
    H2 = new double[3][3][3];
    int i, j, k, l;
    for (i = 0; i < 3; ++i) {
      for (j = 0; j < 3; ++j) {
        for (k = 0; k < 3; ++k) {
          H2[i][j][k] = 0;
          for (l = 0; l < 3; ++l) {
            T[i][j][k][l] = 0;
          }
        }
      }
    }  
    double afGapCost = subMatrix.getAffineGapCost();
    H2[0][1][HZ] = afGapCost;   /* two consecutive horizontal moves */
    H2[1][0][VT] = afGapCost;   /* two consecutive vertical moves */
    T[0][0][0][1] = afGapCost;  /* (-,-) -> (-,x) */  // quasi-natural
    T[0][0][1][0] = afGapCost;  /* (-,-) -> (x,-) */  // quasi-natural
    T[0][1][1][0] = afGapCost;  /* (-,x) -> (x,-) */
    T[1][0][0][1] = afGapCost;  /* (x,-) -> (-,x) */
    T[1][1][0][1] = afGapCost;  /* (x,x) -> (-,x) */
    T[1][1][1][0] = afGapCost;  /* (x,x) -> (x,-) */
  }
    
  public String getSequence(String key) {
    return mapIDtoSeq.get(key).trim();
  }
  
  public int[] getEdgeIDTable() {
    return eIDTab;
  }

  public MSAState initial() {
    // create the initial state                                                                                                                                              
    if (initial == null) {
      initial = new MSAState(numSeqs);
    }
    initial.e = (byte)numOps;
    initial.h = heuristic.getInitH();
    return initial;
  }
  
  public MSANode initialNode() {
    MSAState state = initial();
    MSANode node = new MSANode();
    double h = heuristic.getInitH();
    node.g = 0;
    node.f = h;
    node.packed = pack(state);
    node.e = state.e;
    return node;
  }

  public int expand(MSAState state, MSANode parent, MSANode children[]) {
    return opForward(state, parent, children);
  }
  
  private int opForward(MSAState state, MSANode parent, MSANode children[]) {    
    
    /*
     * Init some essential variables
     */
    int offEdge = 0;            // have we gone off an edge
    boolean reuseCost = true;   // can we reuse the previous cost
    double prevCost = 0;        // the previous cost
        
    /*
     * Compute deltas for incoming edges
     */
    int delta0[] = new int[numSeqs];
    for (int j = 0; j < numSeqs; j++) {
      delta0[j] = (state.e & (1 << j)) > 0 ? 1 : 0;   
    }
    
    /* 
     * Special handling for not penalizing terminal gaps
     */
    if (!penTermGaps) {
      for (int j = 0; j < numSeqs; j++)
        if (state.pos[j] == 0 || state.pos[j] == seqLen[j]) { 
          delta0[j] = 2;
        }
    } 
    
    // Compute holders for the affine costs
    double Tptr[][][] = new double[numSeqs*(numSeqs - 1)/2][3][3];
    for (int i = 1, k = 0; i < numSeqs; i++) {
      for (int j = 0; j < i; j++, k++) {
        Tptr[k][0] = T[delta0[i]][delta0[j]][0];
        Tptr[k][1] = T[delta0[i]][delta0[j]][1];
        Tptr[k][2] = T[delta0[i]][delta0[j]][2];
      }
    }
  
    // Loop over the operators and obtain the gray code
    int generated = 0;
    for (int op=0; op<gcBit.length; op++) {
      int gc = gcBit[op];
      boolean gcInc = gcBitInc[op];
             
      /*
       * increment the seq pos and check offEdge
       */
      if (gcInc) {
        if (state.pos[gc] == seqLen[gc]) offEdge++;
        state.pos[gc]++;
      }
      else {
        state.pos[gc]--;
        if (state.pos[gc] == seqLen[gc]) --offEdge;
      }
      if (offEdge > 0) {
        reuseCost = false; 
        continue; 
      }
        
      // Compute the cost
      int opDim = oprDim[op];
      double cost = 0;
      if (!reuseCost || !gcInc) {
        reuseCost = true;
        char seq[][] = oprSeq[op];
        int pos[] = oprSeqPos[op];
        for (int i = 1; i < opDim; ++i) {
          int p = state.pos[pos[i]]; 
          char ch = seq[i][p];
          for (int j = 0; j < i; j++) {
            int q = state.pos[pos[j]];  
            cost += D[ch][seq[j][q]];
          }
        }
      }
      else {
        char seq[][] = oprSeqPw[op];
        int pos[] = oprSeqPwPos[op];
        char ch = seqs[gc][state.pos[gc]]; 
        for (int j = 0; j < opDim - 1; j++) {
          int p = state.pos[pos[j]]; 
          cost += D[ch][seq[j][p]];    
        }
        cost += prevCost;
      }
      prevCost = cost;
      
      
      // Add affine gap cost 
      int eID = eIDTab[op];
      int delta1[] = new int[numSeqs];
      for (int j = 0; j < numSeqs; j++)
        delta1[j] = (eID & (1 << j)) > 0 ? 1 : 0;
      for (int i = 1, k = 0; i < numSeqs; i++) {
        for (int j = 0; j < i; j++, k++) {
          cost += Tptr[k][delta1[i]][delta1[j]];
        }
      }
      
      // Special handling for not penalizing terminal gaps
      if (!penTermGaps) {
        for (int j = 0; j < numSeqs; j++)
          if (state.pos[j] == 0 || state.pos[j] == seqLen[j]) { 
            delta1[j] = 2;
          }
      } 
      
      // Generate the child nodes
      // TODO: some of this code should be part of the search algorithm?
      double h = 0;
      MSANode child = children[generated];
      try {
        h = heuristic.getH(state, delta1);
        child.parent = parent;
        child.g = parent.g + cost + oprGapCost[op];
        child.f = child.g + h;
        child.packed = pack(state);
        child.e = (byte)eID;
        generated++;
      } catch (Exception e) {
        e.printStackTrace();
      }
      
      // Check for a goal
      if (h == 0) {
        int i = 0;
        for (; i < numSeqs; i++) {
            if (state.pos[i] != seqLen[i]) {                      
                break;
            }
        }
        if (i == numSeqs) {
          setGoal(child);
          return generated;
        }
      }
    }
    
    return generated;
  }
  
  private synchronized void setGoal(MSANode node) {
    if (goal == null || goal.g > node.g) {
      goal = node;
    }
  }

  public MSAState copy(MSAState state) {
    MSAState copy = new MSAState(numSeqs);
    System.arraycopy(state.pos, 0, copy.pos, 0, state.pos.length);
    copy.h = state.h;
    copy.e = state.e;
    return copy;
  }

  public long pack(MSAState state) {
    // this assumes less than 6 sequences
    long word = 0;
    for (int i = 0; i < numSeqs; i++) {
      word = (word << 10) | state.pos[i];
    }
    return word;
  }

  public MSAState unpack(long packed, byte e) {
    MSAState state = new MSAState(numSeqs);
    unpack(packed, e, state);
    return state;
  }

  public void unpack(long packed, byte e, MSAState state) {
    for (int i = numSeqs - 1; i >= 0; i--) {
      int t = (int) packed & 0x3FF;
      packed >>= 10;
      state.pos[i] = t;
    }
    state.e = e;
  }
  
  public void initAbstraction(int abs, int absMod) {
    projection = new MSALatticeProjection(abs, absMod, 
        longestSeqLength, longestSeqLength2, longestSeqIndex, longestSeqIndex2);
  }

  public int getProjectionInt(MSAState state) {
    return (int)projection.project(state);
  }

  public int getProjectionSizeInt() {
    return (int)projection.size();
  }
  
  public long getProjection(MSAState state) {
    return projection.project(state);
  }

  public long getProjectionSize() {
    return projection.size();
  }

  public int getRootOp() {
    return numOps+1;
  }
  
  public int getNumOps() {
    return numOps;
  }
  
  public String seqToString() {
    StringBuffer sb = new StringBuffer();
    for (String key : mapIDtoSeq.keySet()) {
      sb.append(key+"\n"+mapIDtoSeq.get(key)+"\n");
    }
    return sb.toString();
  }
  
  public String opsToString() {
    StringBuffer sb = new StringBuffer();
    
    sb.append("\neIDTab gcBit gcBitInc opDim oprGapCost:\n");
    for (int o = 0; o<numOps; o++) {
      int v = eIDTab[o];
      char binstr[] = new char[17];
      binstr[16] = '\0' ;
      for (int i=0; i<16; i++) {
        binstr[15-i] = ((v & 1) == 1) ? '1' : '0' ;
        v = v / 2 ;
      }
      String bits = new String(binstr).trim();
      sb.append(bits);
      sb.append(" "+gcBit[o]);
      sb.append(" "+(gcBitInc[o] ? "1" : "0"));
      sb.append(" "+oprDim[o]);
      sb.append(" "+(int)oprGapCost[o]+"\n");
    }
    
    sb.append("\noprSeq oprSeqPw oprSeqPos oprSeqPwPos:\n");
    for (int o = 0; o<numOps; o++) {
      for (int i=0; i<numSeqs; i++) {
        sb.append(indexSeq(oprSeq[o][i]));
        sb.append(" "+oprSeqPos[o][i]+" ");
        sb.append(indexSeq(oprSeqPw[o][i]));
        sb.append(" "+oprSeqPwPos[o][i]+"\n");
      }
    }
    
    return sb.toString();
  }
  
  private int indexSeq(char[] c) {
    if (c == null) return -1;
    String id = null;
    String seq = new String(c);
    for (String key : mapIDtoSeq.keySet()) {
      if (seq.equals(mapIDtoSeq.get(key))) {
        id = key; break;
      }
    }
    for (int key : mapIndextoID.keySet()) {
      if (id.equals(mapIndextoID.get(key))) {
        return key;
      }
    }
    return -1;
  }
  
  public String gapTableToString() {
    StringBuffer sb = new StringBuffer();
    for (int i=0; i<3; i++) {
      for (int j=0; j<3; j++) {
        for (int k=0; k<3; k++) {
          for (int l=0; l<3; l++) {
            sb.append((int)T[i][j][k][l]+"\n");
          }
        }
      }
    }
    sb.append("-\n");
    for (int i=0; i<3; i++) {
      for (int j=0; j<3; j++) {
        for (int k=0; k<3; k++) {
          sb.append((int)H2[i][j][k]+"\n");
        }
      }
    }
    return sb.toString();
  }
  
  public String alignmentToString() {
    StringBuffer sb = new StringBuffer();
    StringBuffer seqb[] = new StringBuffer[numSeqs];
    for (int i=0; i<seqb.length; i++) {seqb[i] = new StringBuffer();}
    ArrayList<MSANode>path = new ArrayList<MSANode>();
    for (MSANode p = goal; p != null; p = p.parent) {
      path.add(0, p);
    }
    
    MSANode node = path.get(0);
    MSAState parent = unpack(node.packed, node.e);
    for (int p=1; p<path.size(); p++) {
      node = path.get(p);
      MSAState state = unpack(node.packed, parent.e);
      for (int i=0; i<numSeqs; i++) {
        if (parent != null && 
            parent.pos[i] == state.pos[i]) {
          seqb[i].append("-");
        }
        else {
          seqb[i].append(seqs[i][state.pos[i]]);
        }
      }
      parent = state;
    }
    for (int i=0; i<seqb.length; i++) {
      sb.append(seqb[i]+"\n");
    }
    return sb.toString();
  }
  
  public String alignmentToMSFString() {
    StringBuffer sb = new StringBuffer();
    StringBuffer seqb[] = new StringBuffer[numSeqs];
    for (int i=0; i<seqb.length; i++) {seqb[i] = new StringBuffer();}
    ArrayList<MSANode>path = new ArrayList<MSANode>();
    for (MSANode p = goal; p != null; p = p.parent) {
      path.add(0, p);
    }
    
    // print header
    sb.append("//\n\n"); 
    
    MSANode node = path.get(0);
    MSAState parent = unpack(node.packed, node.e);
    for (int p=1; p<path.size(); p++) {
      node = path.get(p);
      MSAState state = unpack(node.packed, parent.e);
      for (int i=0; i<numSeqs; i++) {
        if (parent != null && 
            parent.pos[i] == state.pos[i]) {
          seqb[i].append("-");
        }
        else {
          seqb[i].append(seqs[i][state.pos[i]]);
        }
      }
      parent = state;
    }
    
    int pos = 0;
    while (pos < seqb[0].length()) {
      for (int i=0; i<seqb.length; i++) {
        sb.append(mapIndextoID.get(i)+"\t");
        for (int j=0; j<50; j+=10) {
          if (pos+j >= seqb[i].length()) break;
          sb.append(seqb[i].substring(pos+j, Math.min(seqb[i].length(), pos+j+10)));
          sb.append(" ");
        }
        sb.append("\n");
      }
      sb.append("\n\n");
      pos+=50;
    }
    return sb.toString();
  }
  
  public MSAHeuristic getHeuristic() {
    return this.heuristic;
  }
  
  public int getNumSeqs() {
    return numSeqs;
  }
  
  public int getLongestSeqLength() {
    return longestSeqLength;
  }
  
  public String getLongestSeqID() {
    return mapIndextoID.get(longestSeqIndex);
  }
  
  public SubMatrix getSubMatrix() {
    return subMatrix;
  }
  
}
