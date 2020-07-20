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

import java.io.*;
import java.util.*;

//Handles multiple tabix indexed files / stores matrix [ variant -> samples ]
//Varargs experiment with passing multiple parameters
public class MultiVCFHandler {
  List<String> vcffiles;
  VariantSampleMatrix variantmatrix;
  //List<String> variantkey; //consist of unique keys of variants (alternate alleles are separated); chr_pos_ref_alt
  //problem with this argument
  //public MultiVCFHandler(String... vcfs) throws IOException, VCFReader.FormatException{
  public MultiVCFHandler(List<String> vcfs) throws IOException, VCFReader.FormatException{
    //this.vcflist=new ArrayList<String>();
    //this.variantkey=new ArrayList<String>();
    List<String> vcflist=new ArrayList<String>();
    this.vcffiles=new ArrayList<String>();
    List<String> variantkey=new ArrayList<String>();
    this.variantmatrix=new VariantSampleMatrix();
    if(vcfs.size() > 0){
      //Phase 1: reading in all the VCFs and creating variant keys for all
      //Phase 2: filling up the occurence matrix
      for(int j=0; j < vcfs.size(); j++){
        String v=vcfs.get(j);
        vcflist.add(v);
        //this.variantmatrix.addSample(v);
        //System.out.println(v);
        boolean gowhile=true;
        VCFReader vreader=new VCFReader(v);
        //System.out.println(v);
        while(gowhile){
          try{
            VCFReader.Variant var=vreader.getNextVariant();
            for(String vkey:var.getVariantKeys()){
              if ( ! this.variantmatrix.getVariations().contains(vkey)) {
                variantkey.add(vkey);
              }
            }
          }
          catch(NullPointerException ne){
            gowhile=false;
          }
        }
        vreader.close();
      }
      //initialize the variant matrix and fill it using the variants in the files
      this.vcffiles=vcflist;
      this.variantmatrix=new VariantSampleMatrix(variantkey, vcflist);
      for(int v=0; v < vcflist.size(); v++){ //column
        TabixReader tr=new TabixReader(vcflist.get(v)+".gz");
        String s;
        for(int i=0; i < this.variantmatrix.getNumberOfVariants(); i++){
          String[] skey=this.variantmatrix.getVariantKey(i).split("#"); //chr_pos_ref_alt
          try{
            TabixReader.Iterator iter=tr.query(skey[0]+":"+skey[1]+"-"+skey[1]);
            while(iter!=null && (s=iter.next()) != null ){
              VCFReader.Variant gvar=new VCFReader.Variant(s);
              for (int gi=0; gi < gvar.getNumberOfAlternateAlleles() ; gi++ ) {
                if(skey[0].equals( gvar.getChromosome() ) && Integer.parseInt(skey[1])==gvar.getPosition() && skey[2].equals(gvar.getReferenceAllele()) && skey[3].equals( gvar.getAlternateAllele(gi) ) ){
                  //do stuff
                  this.variantmatrix.setElement(i,v,1);
                }
              }
            }
          }
          catch(ArrayIndexOutOfBoundsException ae){
            System.out.println("Ignored variant: "+this.variantmatrix.getVariantKey(i));
          }
        }
      }
      System.out.println("Completed processing : "+vcfs.size());
    }
    System.out.println("Matrix dimentions: Rows, Columns - "+Integer.toString(this.variantmatrix.getNumberOfRows())+" , "+Integer.toString(this.variantmatrix.getNumberOfColumns()));
  }
  int getNumberOfVCFFiles(){
    return this.variantmatrix.getNumberOfSamples();
  }
  int getNumberOfSamples(){
    return this.getNumberOfVCFFiles();
  }
  int getNumberOfVariants(){
    return this.variantmatrix.getVariations().size();
  }
  void writeAnnovarInputFormat(String fi) throws IOException, IndexOutOfBoundsException{
    FileWriter fw=new FileWriter(fi);
    fw.write("#chrom\tstartpos\tendpos\tref\talt\n");
    for(int i=0; i < this.getNumberOfVariants(); i++){
      String[] svkey=this.variantmatrix.getVariantKey(i).split("#");
      int end=0;
      end=Integer.parseInt(svkey[1]);
      String endpos=Integer.toString(end);
      if (svkey[2].length()>1) {
        endpos= Integer.toString(end+svkey[2].length()-1);
      }
      fw.write(svkey[0]+"\t"+svkey[1]+"\t"+endpos+"\t"+svkey[2]+"\t"+svkey[3]+"\n");
    }
    fw.close();
  }
  VariantSampleMatrix getVariantSampleMatrix(){
    return this.variantmatrix;
  }
  List<String> getVCFFileList(){
    return this.variantmatrix.getSamples();
  }
  List<String> getUniqueVariantKeys(){
    return this.variantmatrix.getVariations();
  }
  private int indexOfVCFFile(String vfile){
    return this.variantmatrix.getSamples().indexOf(vfile);
  }
  private int indexOfVariantKey(String vkey){
    return this.variantmatrix.getVariations().indexOf(vkey);
  }
  int getNumberOfVariantsPerSample(String vid){
    return this.variantmatrix.getSampleMutationCount(vid);
  }
  int getVariantRecurrence(String vkey){
    return this.variantmatrix.getVariantRecurrence(vkey);
  }
  int getVariantRecurrence(VCFReader.Variant v){
    return this.variantmatrix.getVariantRecurrence(v);
  }
  List<String> getGenomicIntervalVariantKeys(GenomicInterval gi){
    List<String> rvkeys=new ArrayList<String>();
    for(String vkey:this.variantmatrix.getVariations()){
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
  int getTotalVariationCount(GTFReader gr, String sizesfi) throws IOException, GTFReader.FormatException, GTFParseException, GenomicInterval.GenomicIntervalException{
    int totvar=0;
    GenomicIntervalList gil=new GenomicIntervalList(gr);
    GenomicIntervalList ugil=gil.getNonOverlappingGenomicIntervalList(sizesfi);
    for (int i=0; i < ugil.getNumberOfGenomicIntervals() ; i++ ) {
        totvar+= this.getNumberOfVariations(ugil.getGenomicInterval(i));
    }
    return totvar;
  }
  void writeVariantMatrix(String fi) throws IOException{
    FileWriter fw=new FileWriter(fi);
    fw.write("\t");
    System.out.println("Number of samples  : "+Integer.toString(this.variantmatrix.getNumberOfSamples()));
    System.out.println("Number of variants : " + Integer.toString(this.variantmatrix.getNumberOfVariants()));
    for(int i=0; i < this.variantmatrix.getNumberOfColumns(); i++)
    {
      fw.write("\t"+this.variantmatrix.getSample(i));
    }
    fw.write("\n");
    for(int r=0; r< this.variantmatrix.getNumberOfRows(); r++){
      fw.write(this.variantmatrix.getVariantKey(r));
      for(int c=0; c < this.variantmatrix.getNumberOfColumns(); c++){
        fw.write("\t"+Integer.toString(this.variantmatrix.getElement(r,c)));
      }
      fw.write("\n");
    }
  }
}
