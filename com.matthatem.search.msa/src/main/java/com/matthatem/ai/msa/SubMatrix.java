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
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;

/**
 * The substitution matrix class class.
 * 
 * @author Matthew Hatem
 */
public class SubMatrix {
  
  private static final char DASH = '-';
  private static final int MAX_COST_ENTRY = 1000;
  private static final double DEFAULT_F_BOUND = Integer.MAX_VALUE/4;
  private static final double DEFAULT_LINEAR_GAP_COST = 2;
  private static final double DEFAULT_AFFINE_GAP_COST = 8;
  
  private double lrGapCost = DEFAULT_F_BOUND;
  private double afGapCost = DEFAULT_AFFINE_GAP_COST;
  private double tmGapCost;
  private boolean penTermGap;
  
  public double D[][] = new double[256][256];
  char[] ch1 = new char[MAX_COST_ENTRY];
  char[] ch2  = new char[MAX_COST_ENTRY];
  
  public SubMatrix(String path, boolean penTermGap) {
    this.penTermGap = penTermGap;
    File file = new File(path);
    InputStream stream;
    try {
      stream = new FileInputStream(file);
      BufferedReader reader = new BufferedReader(
          new InputStreamReader(stream));
      initCostMatrix(reader);
      stream.close();
    } catch (FileNotFoundException e) {
      e.printStackTrace();
    } catch (IOException e) {
      e.printStackTrace();
    }
  }

  private void initCostMatrix(BufferedReader reader) throws IOException {
    double[] cost = new double[MAX_COST_ENTRY];
    double DMin = DEFAULT_F_BOUND;
    int entry = 0;
    
    /*
     * Read a big 3D array from file
     */
    String nextLine = null;
    while ((nextLine = reader.readLine()) != null) {
      String[] tokens = nextLine.split(" ");
      switch (tokens[0].toUpperCase().charAt(0)) {    
      case DASH: /* DASH ('-') is a reserved character for affine gap cost */
        afGapCost = Double.parseDouble(tokens[1]);
        break;
      default:
        ch1[entry] = tokens[0].toUpperCase().charAt(0);
        ch2[entry] = tokens[1].toUpperCase().charAt(0);
        cost[entry] = Double.parseDouble(tokens[2]);
        if (ch2[entry] == DASH) {
          if (lrGapCost == DEFAULT_F_BOUND) {
            lrGapCost = cost[entry];
          }
          else if (cost[entry] != lrGapCost) {
            System.err.println("gap cost must be linear!\n");
            System.exit(1);
          }
        }
        if (DMin > cost[entry])
          DMin = cost[entry];
        ++entry;      
        break;
      }
    }
    
    /*
     * Build the matrix and default values
     */
    if (lrGapCost == DEFAULT_F_BOUND)
      lrGapCost = DEFAULT_LINEAR_GAP_COST;
    for (int i = 0; i < 256; ++i) 
      for (int j = 0; j < 256; ++j)
        D[i][j] = 255;
    double d;
    for (int i = 0; i < entry; ++i) {
      if (ch2[entry] != DASH)
        d = cost[i] - 2 * DMin;     
      else 
        d = cost[i] - DMin;
      
      D[toupper(ch1[i])][toupper(ch2[i])] = d;
      D[toupper(ch2[i])][toupper(ch1[i])] = d;
      D[tolower(ch1[i])][tolower(ch2[i])] = d;
      D[tolower(ch2[i])][tolower(ch1[i])] = d;
      
      D[toupper(ch1[i])][tolower(ch2[i])] = d;
      D[toupper(ch2[i])][tolower(ch1[i])] = d;
      D[tolower(ch1[i])][toupper(ch2[i])] = d;
      D[tolower(ch2[i])][toupper(ch1[i])] = d;

    }
    lrGapCost -= DMin;
    afGapCost -= DMin;
    if (penTermGap) {
      tmGapCost = afGapCost;
    }
    else {
      tmGapCost = 0;
    }
  }
  
  public char toupper(char c) {
    return Character.toUpperCase(c);
  }
  
  public char tolower(char c) {
    return Character.toLowerCase(c);
  }
  
  public double[][] getTable() {
    return D;
  }
  
  public double getLinearGapCost()  {
    return lrGapCost;
  }
  
  public double getAffineGapCost() {
    return afGapCost;
  }
  
  public double getTerminalGapCost() {
    return tmGapCost;
  }

  public String toString() {
    StringBuffer sb = new StringBuffer();
    sb.append("        "+ ch2[0]);
    for (int i=1; i<27; i++) {
      sb.append("    "+ ch2[i]);
    }
    sb.append("\n");
    for (int i=0; i<27; i++) {
      sb.append("   "+ch2[i]);
      for (int j=0; j<27; j++) {
        sb.append("   "+(int)D[ch2[i]][ch2[j]]);
      }
      sb.append("\n");
    }
    return sb.toString();
  }
    
}