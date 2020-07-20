/*
   The MIT License

   Copyright (c) 2020 Sanket Desai.

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
 * Email    : desai.sanket12@gmail.com
*/

package jabix;
import java.io.*;
import java.util.*;

public class SAMHeaderRecord{
  Map<String, String> header_map;
  String header_type;
  public SAMHeaderRecord(){}
  public SAMHeaderRecord(String se) throws SAMHeaderRecordException, ArrayIndexOutOfBoundsException {
    this.header_map=new HashMap<String, String>();
    this.header_type="";
    if( se.startsWith("@") ){
      String[] splitse=se.split("\t");
      if(splitse.length < 1){
        throw new SAMHeaderRecordException("Header line does not contain data!! Please check!!");
      }
      else{
        this.header_type=splitse[0].substring(1,splitse[0].length());
        if (splitse.length > 1) {
          for(String sps : Arrays.copyOfRange(splitse, 1, splitse.length ) ){ //Starting from 1 to end of elements of the array
            if(sps.contains(":")){
              String[] stse=sps.split(":");
              if(stse.length>1){
                this.header_map.put(stse[0], stse[1]);
              }
            }
            else{
              this.header_map.put(sps, "");
            }
          }
        }
      }
    }
    else{
      throw new SAMHeaderRecordException("Line does not start with '@', hence not a SAM header record!!");
    }
    /*this.salignmentrecord_=se.split("\t");
    if (this.salignmentrecord_.length < 1 ) {
      throw new SAMHeaderRecordException("Number of fields are less in the SAM Header Record. Please refer SAM format documentation!!");
    }
    */
  }
  boolean isHeaderType(String ht){
    return this.header_type.equals(ht);
  }
  String getHeaderType(){
    return this.header_type;
  }

  String getHeaderFeatureValue(String key){
    return this.header_map.get(key);
  }

  public String toString(){
    String header_rec="@"+this.header_type;
    for(Map.Entry <String, String> me:this.header_map.entrySet()){
      header_rec=header_rec+"\t"+me.getKey()+":"+me.getValue();
    }
    return header_rec;
  }
  public class SAMHeaderRecordException extends Exception {
    public SAMHeaderRecordException(){ super("SAM Record format exception. Please refer http://samtools.github.io/hts-specs/SAMv1.pdf ."); }
    public SAMHeaderRecordException(String er){ super(er); }
  }
}
