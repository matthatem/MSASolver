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

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Date;
import java.util.List;
import java.util.Properties;

/**
 * The search result class.
 * 
 * @author Matthew Hatem
 *
 * @param <T> the state type
 */
public class SearchResult<T> {
  
  private List<T> path;
  private Properties properties = new Properties();
  private long startWall;
  private long endWall;
  private long expanded;
  private long generated;
  private int initH;
  private double cost;
  private String algorithm;
   
  /**
   * The constructor.
   */
  public SearchResult() {
  }
  
  /**
   * The constructor.
   * 
   * @param path the solution path
   * @param expanded the number of nodes expanded
   * @param generated the number of nodes generated
   */
  public SearchResult(List<T> path, long expanded, long generated) {
    this.path = path;
    this.expanded = expanded;
    this.generated = generated;
  }
  
  public void setProperty(String key, String value) {
    properties.setProperty(key, value);
  }
  
  public void setProperty(String key, int value) {
    properties.setProperty(key, Integer.toString(value));
  }
    
  /**
   * Sets the solution path.
   * 
   * @param path the solution path
   */
  public void setPath(List<T> path) {
    this.path = path;
  }
  
  /**
   * Sets the number of nodes expanded.
   * 
   * @param expanded the number of nodes expanded
   */
  public void setExpanded(long expanded) {
    this.expanded = expanded;
  }
  
  public double getExpanded() {
    return this.expanded;
  }
  
  /**
   * Sets the number of nods generated.
   * 
   * @param generated the number of nodes generated
   */
  public void setGenerated(long generated) {
    this.generated = generated;
  }
    
  public double getGenerated () {
    return this.generated;
  }
  
  /**
   * Sets the initial heuristic value.
   * 
   * @param initH the initial heuristic
   */
  public void setInitialH(int initH) {
    this.initH = initH;
  }
  
  /**
   * Returns the solution path.
   * 
   * @return the solution path
   */
  public List<T> getPath() {
    return path;
  }
  
  /**
   * Sets the start time in milliseconds.
   * 
   * @param start the start time
   */
  public void setStartTime(long start) {
    this.startWall = start;
  }
  
  /**
   * Sets the end time in milliseconds.
   * 
   * @param end the end time
   */
  public void setEndTime(long end) {
    this.endWall = end;
  }

  /**
   * Returns the machine Id.
   * 
   * @return the machine id
   */
  public String getMachineId() {
    String uname = "unknown";
    try {
      String switches[] = new String[] {"n", "s", "r", "m"};
      String tokens[] = new String[4];
      for (int i=0; i<switches.length; i++) {
        Process p = Runtime.getRuntime().exec("uname -"+switches[i]);
        p.waitFor();
        BufferedReader reader = new BufferedReader(
            new InputStreamReader(p.getInputStream()));
        tokens[i] = reader.readLine();
      }
      uname = tokens[0]+"-"+tokens[1]+"-"+tokens[2]+"-"+tokens[3];
    } catch (IOException e) {
      e.printStackTrace();
    } catch (InterruptedException e) {
      e.printStackTrace();
    }
    
    return uname;
  }
  
  /**
   * Returns a string representation of the results.
   * 
   * @return a string representation
   */
  public String toString() {
    StringBuffer sb = new StringBuffer();
    sb.append("#start data file format 4\n");
    sb.append("#pair\t\"machine id\"\t\""+getMachineId()+"\"\n");
    sb.append("#pair\t\"algorithm\"\t\""+algorithm+"\"\n");
    sb.append("#pair\t\"total wall time\"\t\""+((endWall-startWall)/1000.0)+"\"\n");
    sb.append("#pair\t\"total nodes expanded\"\t\""+expanded+"\"\n");
    sb.append("#pair\t\"total nodes generated\"\t\""+generated+"\"\n");
    sb.append("#pair\t\"solution length\"\t\""+path.size()+"\"\n");
    sb.append("#pair\t\"solution cost\"\t\""+cost+"\"\n");
    sb.append("#pair\t\"initial heuristic\"\t\""+initH+"\"\n");
    sb.append("#pair\t\"wall start date\"\t\""+new Date(startWall)+"\"\n");
    sb.append("#pair\t\"wall start time\"\t\""+startWall+"\"\n");
    sb.append("#pair\t\"wall finish date\"\t\""+new Date(endWall)+"\"\n");
    sb.append("#pair\t\"wall finish time\"\t\""+endWall+"\"\n");
    
    // optional properties
    for (Object key : properties.keySet()) {
      String value = properties.getProperty((String)key);
      sb.append("#pair\t\""+key+"\"\t\""+value+"\"\n");
    }
    
    sb.append("#end data file format 4\n");
    
    return sb.toString();
  }
  
  public double getCost() {
    return cost;
  }

  public void setCost(double cost) {
    this.cost = cost;
  }  
  
  public void setAlgorithm(String algorithm) {
	this.algorithm = algorithm;
  }

}
