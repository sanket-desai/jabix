/*
   The MIT License

   Copyright (c) 2018 Sanket Desai.

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

/*
 * Developer: Sanket Desai
 * Email    : desai.sanket12@outlook.com
*/

package jabix;

import java.util.*;
import java.io.*;
import java.io.FileWriter; //for write method
public class VariantSampleMatrix extends IntegerMatrix {
  List<String> variations;
  List<String> samples;
  public VariantSampleMatrix(){
    this.variations=new ArrayList<String>();
    this.samples=new ArrayList<String>();
    this.data=new int[0][0];
  }
  public VariantSampleMatrix(int rw, int cl){
    super(rw, cl);
    this.variations=new ArrayList<String>();
    this.samples=new ArrayList<String>();
  }
  public VariantSampleMatrix(List<String> vari, List<String> samp){
    super(vari.size(),samp.size());
    this.variations=vari;
    this.samples=samp;
  }
  public VariantSampleMatrix(String fname) throws IOException, FileNotFoundException {
    this.variations=new ArrayList<String>();
    this.samples=new ArrayList<String>();
    String sep="\t";
    //File containing binary matrix of variants->samples and cell containing integer (mutation presence 1 / 0 abscence)
    BufferedReader br=new BufferedReader( new FileReader(fname) );
    String line=br.readLine();
    String[] sline=line.split(sep);
    //if ( ! sline[0].equals("") ) {
    if (! "".equals(sline[0])){
      System.out.println("Format Error: First cell of the matrix should be blank!! Please check input file!!");
      System.exit(0);
    }
    for(int i=1; i < sline.length; i++){
      this.samples.add(sline[i]);
    }
    List<String> temprows=new ArrayList<String>();
    line=br.readLine();
    temprows.add(line);
    while(line!=null){
      temprows.add(line);
      line=br.readLine();
    }
    this.setNumberOfRows(temprows.size());
    this.setNumberOfColumns(this.samples.size());
    this.data=new int[this.rows][this.cols];
    for(int r=0; r < temprows.size(); r++){
      String[] srow=temprows.get(r).split(sep);
      if (srow.length-1 != this.samples.size()) {
        System.out.println("Format Error: Encountered row having less columns! Please check input file!!");
        System.exit(0);
      }
      this.variations.add(srow[0]);
      for (int c=1; c < srow.length ; c++ ) {
        //srow[x] can exist in two formats; a) 0 or 1 b) 0|1 (separated by '|')
        if (srow[c].contains("|")) {
          String[] gtcell=srow[c].split("|");
          if (gtcell[0]=="1" || gtcell[1]=="1") {
            //this.setElement(r, c-1, Integer.parseInt(srow[c]));
            this.setElement(r, c-1, 1);
          }
          else{
            this.setElement(r,c-1, 0);
          }
        }
        else if (srow[c].contains(",")) {
          String[] gtcell=srow[c].split(",");
          if (gtcell[0]=="1" || gtcell[1]=="1") {
            //this.setElement(r, c-1, Integer.parseInt(srow[c]));
            this.setElement(r, c-1, 1);
          }
          else{
            this.setElement(r,c-1, 0);
          }
        }
        else{
          this.setElement(r, c-1, Integer.parseInt(srow[c]));
        }
      }
    }
  }
  public void addVariation(String v){
    if (!this.variations.contains(v)) {
      this.variations.add(v);
    }
  }
  public void addSample(String s){
    if (!this.samples.contains(s)) {
      this.samples.add(s);
    }
  }
  public void setVariations(List<String> vars){
    this.variations=vars;
  }
  public void setSamples(List<String> sam){
    this.samples=sam;
  }
  public int getNumberOfSamples(){
    return this.samples.size();
  }
  public int getNumberOfVariants(){
    return this.variations.size();
  }
  public List<String> getVariations(){
    return this.variations;
  }
  public String getVariantKey(int i){
    return this.getVariations().get(i);
  }
  public List<String> getSamples(){
    return this.samples;
  }
  public String getSample(int i){
    return this.getSamples().get(i);
  }
  public int getSampleMutationCount(String samplevcf){
    return this.getColumnSum(this.samples.indexOf(samplevcf));
  }
  public int getVariantRecurrence(String vkey){
    return this.getRowSum(this.variations.indexOf(vkey));
  }
  public int getVariantRecurrence(VCFReader.Variant v){
    int reccur=0;
    List<String> vkeys=v.getVariantKeys();
    for(int i=0; i < vkeys.size(); i++){
      int arec=this.getRowSum(this.samples.indexOf(vkeys.get(i)));
      if(arec > reccur)
      {
        reccur=arec;
      }
    }
    return reccur;
  }
  List<String> getGenomicIntervalVariantKeys(GenomicInterval gi){
    List<String> rvkeys=new ArrayList<String>();
    for(String vkey:this.getVariations()){
      String[] svkey=vkey.split("#");
      if (gi.isWithinRange(svkey[0], Long.parseLong(svkey[1]))) {
        rvkeys.add(vkey);
      }
    }
    return rvkeys;
  }
  //Given genomic interval returns variant count for the same (total/ rowSums addition)
  int getNumberOfVariations(GenomicInterval gi){
    List<String> givarkeys=this.getGenomicIntervalVariantKeys(gi);
    int nvars=0;
    for (String var: givarkeys) {
      nvars+=this.getVariantRecurrence(var);
    }
    return nvars;
  }
  public void writeMatrix(String fn) throws IOException {
    FileWriter fw=new FileWriter(fn);
    for(int i=0; i < this.samples.size(); i++)
    {
      fw.write("\t"+this.samples.get(i));
    }
    fw.write("\n");
    for(int r=0; r< this.variations.size(); r++){
      fw.write(this.variations.get(r));
      for(int c=0; c < this.samples.size(); c++){
        fw.write("\t"+Integer.toString(this.getElement(r,c)));
      }
      fw.write("\n");
    }
  }
  /*
  public static void main(String[] argv){
    try{
      VariantSampleMatrix vsm=new VariantSampleMatrix("/run/media/groot/Vault/ACTREC/DuttLab/dry_lab/data/mutation/databases/1000Genome/chrY.genotype.cleaned.mat");
      String vkey=vsm.getVariantKey(54814);
      System.out.println(vsm.getVariantRecurrence(vkey));
    }
    catch(Exception e){
      e.printStackTrace();
    }
  }
  */
}
