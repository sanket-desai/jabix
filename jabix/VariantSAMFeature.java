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
import net.sf.samtools.SAMUtils;

public class VariantSAMFeature{
  private String samtabfile="";
  private VCFReader.Variant variant;
  //public VariantSAMFeature(){}
  //samfn -> sam file indexed using tabix indexer
  public VariantSAMFeature(VCFReader.Variant var, String samfn) throws GenomicInterval.GenomicIntervalException, IOException {
    this.samtabfile=samfn;
    this.variant=var;
  }
  String getSAMFileName(){
    return this.samtabfile;
  }
  int getVariantDepth() throws GenomicInterval.GenomicIntervalException, IOException {
    int depth=0;
    String s;
    TabixReader tr= new TabixReader(this.samtabfile);
    TabixReader.Iterator iter=tr.query(this.variant.toGenomicInterval());
    while (iter != null && (s = iter.next()) != null){
      depth++;
    }
    return depth;
  }
  List<SAMAlignmentRecord> getSAMAlignmentRecords() throws GenomicInterval.GenomicIntervalException, IOException, SAMAlignmentRecord.SAMAlignmentRecordException {
    String s;
    TabixReader tr= new TabixReader(this.samtabfile);
    TabixReader.Iterator iter=tr.query(this.variant.toGenomicInterval());
    List<SAMAlignmentRecord> sar=new ArrayList<SAMAlignmentRecord>();
    while (iter != null && (s = iter.next()) != null){
      sar.add(new SAMAlignmentRecord(s));
    }
    return sar;
  }
  int getNumberOfReferenceAlleleMappingReads() throws GenomicInterval.GenomicIntervalException, IOException, SAMAlignmentRecord.SAMAlignmentRecordException {
    String s;
    int refreads=0;
    TabixReader tr= new TabixReader(this.samtabfile);
    TabixReader.Iterator iter=tr.query(this.variant.toGenomicInterval());
    while (iter != null && (s = iter.next()) != null){
      SAMAlignmentRecord sar=new SAMAlignmentRecord(s);
      int sindexstart= this.variant.getPosition()-sar.getAlignedPosition(); //where on SAM record
      int sindexend=sindexstart+this.variant.getReferenceAllele().length()  ; //where on SAM record
      String samseq=sar.getReadSequence().substring(sindexstart, sindexend);
      if (this.variant.getReferenceAllele().equals(samseq)) {
        refreads++;
      }
    }
    return refreads;
  }
  List<SAMAlignmentRecord> getReferenceAlleleMappingReads() throws GenomicInterval.GenomicIntervalException, IOException, SAMAlignmentRecord.SAMAlignmentRecordException {
    String s;
    List<SAMAlignmentRecord> refreads=new ArrayList<SAMAlignmentRecord>();
    TabixReader tr= new TabixReader(this.samtabfile);
    TabixReader.Iterator iter=tr.query(this.variant.toGenomicInterval());
    while (iter != null && (s = iter.next()) != null){
      SAMAlignmentRecord sar=new SAMAlignmentRecord(s);
      int sindexstart= this.variant.getPosition()-sar.getAlignedPosition(); //where on SAM record
      int sindexend=sindexstart+this.variant.getReferenceAllele().length()  ; //where on SAM record
      String samseq=sar.getReadSequence().substring(sindexstart, sindexend);
      if (this.variant.getReferenceAllele().equals(samseq)) {
        refreads.add(sar);
      }
    }
    return refreads;
  }
  int getNumberOfAlternateAlleleMappingReads() throws GenomicInterval.GenomicIntervalException, IndexOutOfBoundsException, IOException, SAMAlignmentRecord.SAMAlignmentRecordException {
    String s;
    int altreads=0;
    TabixReader tr= new TabixReader(this.samtabfile);
    TabixReader.Iterator iter=tr.query(this.variant.toGenomicInterval());
    while (iter != null && (s = iter.next()) != null){
      SAMAlignmentRecord sar=new SAMAlignmentRecord(s);
      int sindexstart= this.variant.getPosition()-sar.getAlignedPosition(); //where on SAM record
      int sindexend=sindexstart+this.variant.getAlternateAllele(0).length()  ; //where on SAM record
      String samseq=sar.getReadSequence().substring(sindexstart, sindexend);
      if ( this.variant.getAlternateAllele(0).equals(samseq)) {
        altreads++;
      }
    }
    return altreads;
  }
  List<SAMAlignmentRecord> getAlternateAlleleMappingReads() throws GenomicInterval.GenomicIntervalException, IndexOutOfBoundsException, IOException, SAMAlignmentRecord.SAMAlignmentRecordException {
    String s;
    List<SAMAlignmentRecord> altreads=new ArrayList<SAMAlignmentRecord>();
    TabixReader tr= new TabixReader(this.samtabfile);
    TabixReader.Iterator iter=tr.query(this.variant.toGenomicInterval());
    while (iter != null && (s = iter.next()) != null){
      SAMAlignmentRecord sar=new SAMAlignmentRecord(s);
      int sindexstart= this.variant.getPosition()-sar.getAlignedPosition(); //where on SAM record
      int sindexend=sindexstart+this.variant.getAlternateAllele(0).length()  ; //where on SAM record
      String samseq=sar.getReadSequence().substring(sindexstart, sindexend);
      if ( this.variant.getAlternateAllele(0).equals(samseq)) {
        altreads.add(sar);
      }
    }
    return altreads;
  }
  //Using fastqToPhred from net.sf.samtools.SAMUtils
  List<Integer> getReferenceAlleleBaseQualities() throws GenomicInterval.GenomicIntervalException, IOException, SAMAlignmentRecord.SAMAlignmentRecordException {
    List<Integer> basequal=new ArrayList<Integer>();
    String s;
    TabixReader tr= new TabixReader(this.samtabfile);
    TabixReader.Iterator iter=tr.query(this.variant.toGenomicInterval());
    while (iter != null && (s = iter.next()) != null){
      SAMAlignmentRecord sar=new SAMAlignmentRecord(s);
      int sindexstart= this.variant.getPosition()-sar.getAlignedPosition(); //where on SAM record
      int sindexend=sindexstart+this.variant.getReferenceAllele().length()  ; //where on SAM record
      String samseq=sar.getReadSequence().substring(sindexstart, sindexend);
      if ( this.variant.getReferenceAllele().equals(samseq)) {
        System.out.println(SAMUtils.fastqToPhred(sar.getQualityString().charAt(sindexstart)));
        basequal.add( SAMUtils.fastqToPhred(sar.getQualityString().charAt(sindexstart)) );
      }
    }
    return basequal;
  }
  int getSumOfReferenceAlleleBaseQualities() throws GenomicInterval.GenomicIntervalException, IOException, SAMAlignmentRecord.SAMAlignmentRecordException {
    int sum=0;
    for ( int temp : this.getReferenceAlleleBaseQualities() ) {
      sum+=temp;
    }
    return sum;
  }
  double getAverageReferenceAlleleBaseQuality() throws GenomicInterval.GenomicIntervalException, IOException, SAMAlignmentRecord.SAMAlignmentRecordException {
    double sum=0, counter=0, avg=0;
    for ( int temp : this.getReferenceAlleleBaseQualities() ) {
      sum+=temp;
      counter++;
    }
    if (counter>0) {
      avg=sum/counter;
    }
    return  avg ;
  }
  List<Integer> getAlternateAlleleBaseQualities() throws GenomicInterval.GenomicIntervalException, IndexOutOfBoundsException, IOException, SAMAlignmentRecord.SAMAlignmentRecordException {
    List<Integer> basequal=new ArrayList<Integer>();
    String s;
    TabixReader tr= new TabixReader(this.samtabfile);
    TabixReader.Iterator iter=tr.query(this.variant.toGenomicInterval());
    while (iter != null && (s = iter.next()) != null){
      SAMAlignmentRecord sar=new SAMAlignmentRecord(s);
      int sindexstart= this.variant.getPosition()-sar.getAlignedPosition(); //where on SAM record
      int sindexend=sindexstart+this.variant.getAlternateAllele(0).length()  ; //where on SAM record
      String samseq=sar.getReadSequence().substring(sindexstart, sindexend);
      if ( this.variant.getAlternateAllele(0).equals(samseq)) {
        basequal.add( SAMUtils.fastqToPhred(sar.getQualityString().charAt(sindexstart)) );
      }
    }
    return basequal;
  }
  int getSumOfAlternateAlleleBaseQualities() throws GenomicInterval.GenomicIntervalException, IndexOutOfBoundsException, IOException, SAMAlignmentRecord.SAMAlignmentRecordException {
    int sum=0;
    for ( int temp : this.getAlternateAlleleBaseQualities() ) {
      sum+=temp;
    }
    return sum;
  }
  double getAverageAlternateAlleleBaseQuality() throws GenomicInterval.GenomicIntervalException, IndexOutOfBoundsException, IOException, SAMAlignmentRecord.SAMAlignmentRecordException {
    double sum=0, counter=0, avg=0;
    for ( int temp : this.getAlternateAlleleBaseQualities() ) {
      sum+=temp;
      counter++;
    }
    if (counter>0) {
      avg=sum/counter;
    }
    return  avg;
  }

  List<Integer> getReferenceAlleleReadMappingQualities() throws GenomicInterval.GenomicIntervalException, IOException, SAMAlignmentRecord.SAMAlignmentRecordException {
    List<Integer> mqsum=new ArrayList<Integer>();
    String s;
    TabixReader tr= new TabixReader(this.samtabfile);
    TabixReader.Iterator iter=tr.query(this.variant.toGenomicInterval());
    while (iter != null && (s = iter.next()) != null){
      SAMAlignmentRecord sar=new SAMAlignmentRecord(s);
      int sindexstart= this.variant.getPosition()-sar.getAlignedPosition(); //where on SAM record
      int sindexend=sindexstart+this.variant.getReferenceAllele().length()  ; //where on SAM record
      String samseq=sar.getReadSequence().substring(sindexstart, sindexend);
      if ( this.variant.getReferenceAllele().equals(samseq)) {
        mqsum.add( sar.getMappingQuality() );
      }
    }
    return mqsum;
  }
  List<Integer> getAlternateAlleleReadMappingQualities() throws GenomicInterval.GenomicIntervalException, IndexOutOfBoundsException, IOException, SAMAlignmentRecord.SAMAlignmentRecordException {
    List<Integer> mqsum=new ArrayList<Integer>();
    String s;
    TabixReader tr= new TabixReader(this.samtabfile);
    TabixReader.Iterator iter=tr.query(this.variant.toGenomicInterval());
    while (iter != null && (s = iter.next()) != null){
      SAMAlignmentRecord sar=new SAMAlignmentRecord(s);
      int sindexstart= this.variant.getPosition()-sar.getAlignedPosition(); //where on SAM record
      int sindexend=sindexstart+this.variant.getAlternateAllele(0).length()  ; //where on SAM record
      String samseq=sar.getReadSequence().substring(sindexstart, sindexend);
      if (  this.variant.getAlternateAllele(0).equals(samseq)) {
        mqsum.add(sar.getMappingQuality());
      }
    }
    return mqsum;
  }
  double getAverageReferenceAlleleMappingQuality() throws GenomicInterval.GenomicIntervalException, IOException, SAMAlignmentRecord.SAMAlignmentRecordException {
    double avg=0, counter=0, total=0;
    for (int mq:this.getReferenceAlleleReadMappingQualities()) {
      total+=mq;
      counter++;
    }
    if (counter>0) {
      avg=total/counter;
    }
    return avg;
  }
  double getAverageAlternateAlleleMappingQuality() throws GenomicInterval.GenomicIntervalException, IndexOutOfBoundsException, IOException, SAMAlignmentRecord.SAMAlignmentRecordException {
    double avg=0, counter=0, total=0;
    for (int mq:this.getAlternateAlleleReadMappingQualities()) {
      total+=mq;
      counter++;
    }
    if (counter>0) {
      avg=total/counter;
    }
    return avg;
  }
  //Returns minimum tail distance of mutations from either of the read ends
  List<Integer> getReferenceAlleleReadTailDistances() throws GenomicInterval.GenomicIntervalException, IndexOutOfBoundsException, IOException, SAMAlignmentRecord.SAMAlignmentRecordException {
    List<Integer> taildist=new ArrayList<Integer>();
    String s;
    TabixReader tr= new TabixReader(this.samtabfile);
    TabixReader.Iterator iter=tr.query(this.variant.toGenomicInterval());
    while (iter != null && (s = iter.next()) != null){
      SAMAlignmentRecord sar=new SAMAlignmentRecord(s);
      int sindexstart= this.variant.getPosition()-sar.getAlignedPosition(); //where on SAM record
      int sindexend=sindexstart+this.variant.getReferenceAllele().length()  ; //where on SAM record
      String samseq=sar.getReadSequence().substring(sindexstart, sindexend);
      if ( this.variant.getReferenceAllele().equals(samseq)) {
        int tail1=this.variant.getPosition()-sar.getAlignedPosition()+1; //Right side of a read
        int tail2=sar.getAlignedPosition()+sar.getReadSequenceLength()-this.variant.getPosition();
        if (tail1 < tail2) {
          taildist.add(tail1);
        }
        else if (tail2 < tail1) {
          taildist.add(tail2);
        }
        else{
          taildist.add(tail1);
        }
      }
    }
    return taildist;
  }
  double getAverageReferenceAlleleReadTailDistance() throws GenomicInterval.GenomicIntervalException, IndexOutOfBoundsException, IOException, SAMAlignmentRecord.SAMAlignmentRecordException {
    double avg=0, total=0, counter=0;
    for (int td:this.getReferenceAlleleReadTailDistances()) {
      total+=td;
      counter++;
    }
    if (counter>0) {
      avg=total/counter;
    }
    return avg;
  }
  List<Integer> getAlternateAlleleReadTailDistances() throws GenomicInterval.GenomicIntervalException, IndexOutOfBoundsException, IOException, SAMAlignmentRecord.SAMAlignmentRecordException {
    List<Integer> taildist=new ArrayList<Integer>();
    String s;
    TabixReader tr= new TabixReader(this.samtabfile);
    TabixReader.Iterator iter=tr.query(this.variant.toGenomicInterval());
    while (iter != null && (s = iter.next()) != null){
      SAMAlignmentRecord sar=new SAMAlignmentRecord(s);
      int sindexstart= this.variant.getPosition()-sar.getAlignedPosition(); //where on SAM record
      int sindexend=sindexstart+this.variant.getAlternateAllele(0).length()  ; //where on SAM record
      String samseq=sar.getReadSequence().substring(sindexstart, sindexend);
      if ( this.variant.getAlternateAllele(0).equals(samseq)) {
        int tail1=this.variant.getPosition()-sar.getAlignedPosition()+1; //Right side of a read
        int tail2=sar.getAlignedPosition()+sar.getReadSequenceLength()-this.variant.getPosition();
        if (tail1 < tail2) {
          taildist.add(tail1);
        }
        else if (tail2 < tail1) {
          taildist.add(tail2);
        }
        else{
          taildist.add(tail1);
        }
      }
    }
    return taildist;
  }
  double getAverageAlternateAlleleReadTailDistance() throws GenomicInterval.GenomicIntervalException, IndexOutOfBoundsException, IOException, SAMAlignmentRecord.SAMAlignmentRecordException {
    double avg=0, total=0, counter=0;
    for (int td:this.getAlternateAlleleReadTailDistances()) {
      total+=td;
      counter++;
    }
    if (counter>0) {
      avg=total/counter;
    }
    return avg;
  }
  public static void main(String[] args){
    //arguments: vcf, sam
    try{
      VCFReader vr=new VCFReader(args[0]);
      VCFReader.Variant v=vr.getNextVariant();
      System.out.println(v.toString());
      VariantSAMFeature vs=new VariantSAMFeature(v , args[1]);
      System.out.println(vs.getAlternateAlleleReadTailDistances().get(0));

      //System.out.println(vs.getVariantDepth());
    }
    catch(Exception e){
      e.printStackTrace();
    }
  }
}
