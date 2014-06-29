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
package com.matthatem.ai.search.applications;

import java.io.FileInputStream;
import java.io.FileNotFoundException;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;

import com.matthatem.ai.msa.MSA;
import com.matthatem.ai.msa.MSA.MSAState;
import com.matthatem.ai.msa.algorithms.Astar;
import com.matthatem.ai.msa.algorithms.SearchAlgorithm;
import com.matthatem.ai.msa.algorithms.SearchResult;


/**
 * The main entry point for the MSA solver.
 * 
 * @author Matthew Hatem
 */
public class MSASolver {
  
  private static String algoString;
  
  public static void main(String[] args) {   
    Options options = createOptions();
    CommandLineParser parser = new PosixParser();
    CommandLine cmd = null;
    try {
      cmd = parser.parse(options, args);
    } catch (ParseException e) {
      e.printStackTrace();
      System.exit(1);
    }
    
    MSA msa = createMSAInstance(cmd);
    SearchAlgorithm algo = createSearchAlgorithm(cmd, msa);
    
    System.gc(); 
    
    long t = System.currentTimeMillis();
    SearchResult<MSAState> result = algo.search();
    long td = System.currentTimeMillis();
        
    result.setAlgorithm(algoString);
    result.setInitialH((int)msa.getHeuristic().getInitH()); 
    result.setStartTime(t);
    result.setEndTime(td);
    
    System.out.println(result); 
    //System.out.println(msa.alignmentToString());
    System.out.println(msa.alignmentToMSFString());
  }
  
  private static Options createOptions() {
    Options options = new Options();
    options.addOption("a", "algorithm", true, "search algorithm");
    options.addOption("i", "instance", true, "path to problem instance");
    options.addOption("w", "weight", true, "weight");
    options.addOption("h", "heuristic", true, "heuristic");
    options.addOption("q", "qgaps", false, "quasi natural gap costs");
    return options;
  }

  private static MSA createMSAInstance(CommandLine cmd) {
    MSA msa = null;
    String path = cmd.getOptionValue("i", null);
    if (path != null && path.length() > 1) {
      try {
        // heuristic
        String h = cmd.getOptionValue("h", "2-fold");
        MSA.HEURISTICS heuristic = MSA.HEURISTICS.H2D;
        if ("divconq".equals(h)) {
          heuristic = MSA.HEURISTICS.HDIVCONQ;
        }
        else if ("divconq_int".equals(h)) {
          heuristic = MSA.HEURISTICS.HDIVCONQ_INT;
        }
        // penalize terminal gaps
        boolean penTermGap = true;
        if (cmd.hasOption("q")) {
          penTermGap = false;
        }
        // weight
        double weight = Double.parseDouble(cmd.getOptionValue("w", "1"));
        msa = new MSA(new FileInputStream(path), heuristic, penTermGap, weight);
      } catch (FileNotFoundException e) {
        e.printStackTrace();
      }
    }
    return msa;
  }
  
  private static SearchAlgorithm createSearchAlgorithm(CommandLine cmd, MSA msa) {
    algoString = cmd.getOptionValue("a");
    SearchAlgorithm algo = null;
    if ("astar".equals(algoString) || "wastar".equals(algoString)) {
      algo = new Astar(msa);
    }
    else {
      fatalError("Unsupported algorithm: "+algoString);
    }    
    return algo;
  }
    
  private static void fatalError(final String message) {
    System.err.println(message);
    System.exit(1);
  }
}
