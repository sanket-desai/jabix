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
import java.util.regex.Pattern;
import java.util.regex.Matcher;


public class GTFReader extends BufferedReader {
  String filename="";
  //MetaData meta;
  long gtfrecordindex=0, gtfrecordstartindex=0;
  public GTFReader(String fname) throws IOException, FileNotFoundException, GTFParseException {
    super(new FileReader(fname));
    this.filename=fname;
    //meta=new MetaData(this);
    //System.out.println(fname);
  }
  boolean hasNextGTFRecord() throws IOException{
    return this.ready();
  }
  GTFRecord getNextGTFRecord() throws IOException, FormatException, GTFParseException {
    String newrec=this.readLine();
    this.gtfrecordindex++;
    return new GTFRecord(newrec);
  }
  /*MetaData getMetaData(){
    return this.meta;
  }*/

  public static class FormatException extends GTFParseException{
    public FormatException(){
      super("GTF file format exception! Please check the file and refer to https://www.gencodegenes.org/gencodeformat.html for more details!!");
    }
    public FormatException(String s){
      super(s);
    }
  }
  /*
  public class MetaData{
    private List<String> metadata;
    int metafindexend=0;
    public MetaData(){}
    public MetaData(GTFReader gr) throws IOException{
      String line;
      while( (line=gr.readLine()) != null ){
        if (line.startsWith("#") ) {
          gr.mark(0);
          metadata.add(line);
        }
        else{
          gr.mark(0);
          break;
        }
      }
    }
    public MetaData(List<String> met){
      this.metadata=met;
    }
    void addMetaDataElement(String s){
      this.metadata.add(s);
    }
    String getMetaDataElement(int i) throws IndexOutOfBoundsException {
      return this.metadata.get(i);
    }
    int getNumberOfMetaDataElements(){
      return this.metadata.size();
    }
  }
  */
  //GTFRecord class
  public static class GTFRecord extends GencodeFeature{
    //private static final Pattern ATTRIBUTE_PATTERN = Pattern
    private final Pattern ATTRIBUTE_PATTERN = Pattern
    .compile("^\\s*(.+)\\s(.+)$");
    public GTFRecord(){
      super();
    }
    public GTFRecord(String line) throws FormatException, GTFParseException {
      /*if (line == null || line.startsWith("#")){
          meta.addMetaDataElement(line);
      }
      */
      String[] fields = line.trim().split("\t");
      if (fields.length < 2) {
        throw new FormatException("File is not as per GENCODE GTF format!!");
      }
      this.seqname = fields[0];
      this.source = fields[1];
      this.featureType = FeatureType.fromString(fields[2]);
      try {
          this.start = Integer.valueOf(fields[3]);
      } catch (NumberFormatException e) {
          throw new FormatException("Invalid integer value for start");
      }
      try {
          this.end = Integer.valueOf(fields[4]);
      }
      catch (NumberFormatException e) {
          throw new FormatException("Invalid integer value for end");
      }
      this.strand = Strand.fromString(fields[6]);
      try {
          this.score = Double.valueOf(fields[5]);
          this.frame = Integer.valueOf(fields[7]);
      }
      catch (NumberFormatException ignored) {}
      this.attributes = new HashMap<String, String>();
      if (fields.length >= 8 && fields[8] != null) {
        String[] attrs = fields[8].split(";");
        for (String variableString : attrs) {
          Matcher m = ATTRIBUTE_PATTERN.matcher(variableString);
          if (m.matches()) {
            String key = m.group(1).trim();
            String val = m.group(2).trim();
            val = val.replaceAll("\"", "");
            if (val.length() > 0) {
              this.attributes.put(key, val);
            }
          }
        }
      }
      this.geneId = this.getAttribute("gene_id");
      this.transcriptId = this.getAttribute("transcript_id");
      this.geneType = this.getAttribute("gene_type");
      this.geneStatus = this.getAttribute("gene_status");
      this.geneName = this.getAttribute("gene_name");
      this.transcriptType = this.getAttribute("transcript_type");
      this.transcriptStatus = this.getAttribute("transcript_status");
      this.transcriptName = this.getAttribute("transcript_name");
      try
      {
        this.level = Integer.valueOf(this.getAttribute("level"));
      }
      catch(NumberFormatException ignored) { }
    }
    int size(){
      return this.getEnd()-this.getStart()+1;
    }
    /*
    String getGeneId(){
      return this.geneId;
    }

    String getFeatureType(){
      return this.featureType;
    }

    String getGeneName(){
      return this.geneName;
    }
    String getTranscriptId(){
      return this.transcriptId;
    }
    String getTranscriptName(){
      return this.transcriptName;
    }
    */
    GenomicInterval asGenomicInterval() throws GenomicInterval.GenomicIntervalException{
      return new GenomicInterval(this.getSeqname(), this.getStart(), this.getEnd());
    }
    boolean hasOverlap(GTFReader.GTFRecord rec) throws GenomicInterval.GenomicIntervalException{
      boolean overlap=false;
      return this.asGenomicInterval().hasOverlap(rec.asGenomicInterval());
    }
  }
}
