/*
   The MIT License

   Copyright (c) 2019 Sanket Desai.

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

public class VCFReader extends BufferedReader {
	String filename="";
	MetaData metadata=new MetaData();
	long variantindex=0, variantstartindex=0;
	public VCFReader(String fname) throws IOException, FormatException, FileNotFoundException {
		super(new FileReader(fname));
		filename=fname;
		metadata=new MetaData(this);
		variantindex=metadata.variantStartLineNumber();
	}
	long getVariantStartLineNumber(){
		return this.variantindex;
	}
	MetaData getMataData(){
		return this.metadata;
	}
	boolean hasNextVariant() throws IOException{
		return this.ready();
	}
	Variant getNextVariant() throws IOException, FormatException {
		String newrec=this.readLine();
		this.variantindex++;
		return new Variant(newrec);
	}

	/*  Classes under VCF Metaclass:
	 * 	1 Header
	 * 	2 InfoField
	 * 	3 FormatField
	 * 	4 FilterField
	 * 	5 Exception
	 * 	6 AlternateAllele (ALT)
	 *	7 Pedigree
	 *	8 AssemblyField
	 *	9 ReferenceField
	 *	10 ContigField
	 *	11 FileFormat
	 *	12 Sample
	 */
	public static class FormatException extends Exception{
	public FormatException(){ super("VCF File format (4.2) exception! Please check the file!!");}
	public FormatException(String e){ super(e);}
	}

	public class InfoField{
		Map<String, String> infomap=new LinkedHashMap<String, String>();
		public InfoField(){}
		public InfoField(String inf) throws FormatException{
			if(inf.startsWith("<ID") && inf.endsWith(">")){
				inf=inf.substring(1,inf.length()-1);
				String[] sinf=inf.split(",");
				for(String f:sinf){
					String[] sf=f.split("=");
					try{
						infomap.put(sf[0].trim(),sf[1].trim());
					}
					catch(Exception e){
						throw new FormatException("Error due to : "+f);
					}
				}
			}
			else if(inf.startsWith("##INFO")){
				inf=inf.substring(8,inf.length()-1);
				String[] sinf=inf.split(",");
				for(String f:sinf){
					String[] sf=f.split("=");
					if(sf.length>2){ throw new FormatException("Error due to : "+f); }
					try{
						infomap.put(sf[0].trim(),sf[1].trim());
					}
					catch(Exception e){
						throw new FormatException("Error due to : "+f);
					}
				}
			}
			else{
				throw new FormatException("Error due to : "+inf);
			}
		}
		String getProperty(String k){
			try{
				return this.infomap.get(k);
			}
			catch (Exception e){
				System.err.println("Key does not exist : "+k);
				return "";
			}
		}
		String infoFieldToString(){
			String infield="##INFO=<ID="+this.getProperty("ID")+",Number="+this.getProperty("Number")+",Type="+this.getProperty("Type")+",Description="+this.getProperty("Description")+">";
			return infield;
		}
	}

	public class FileFormatField{
		String fileformat="";
		public FileFormatField(){}
		public FileFormatField(String f) throws FormatException{
			if(f.startsWith("##fileformat")){
				fileformat=f.substring(13).trim();
			}
			else{
				throw new FormatException("Error due to : "+f);
			}
		}
		String getFormat(){
			return this.fileformat;
		}
	}

	public class SourceField{
		String source="";
		public SourceField(){}
		public SourceField(String f) throws FormatException{
			if(f.startsWith("##source")){
				source=f.substring(9).trim();
			}
			else{
				throw new FormatException("Error due to : "+f);
			}
		}
		String getSource(){
			return this.source;
		}
	}

	public class Header{
		List<String> header=new ArrayList<String>();
		public Header(){}
		public Header(String h) throws FormatException{
			if(h.startsWith("#CHROM")){
				h=h.substring(1).trim();
				this.header=Arrays.asList(h.split("\t"));
			}
			else{
				throw new FormatException("Error due to : "+h);
			}
		}
		List<String> getHeaders(){
			return this.header;
		}
		String getHeader(int i){
			return this.header.get(i);
		}
	}

	public class FilterField{
		Map<String, String> filtermap=new LinkedHashMap<String, String>();
		public FilterField(){}
		public FilterField(String f) throws FormatException{
			if(f.startsWith("##FILTER")){
				f=f.substring(10,f.length()-1);
				String[] sf=f.split(",");
				for(String f1:sf){
					String[] sf1=f1.split("=");
					if(sf1.length>2){ throw new FormatException("Error due to : "+f1);}
					try{
						this.filtermap.put(sf[0].trim(),sf[1].trim());
					}
					catch(Exception e){
						throw new FormatException("Error due to : "+f1);
					}
				}
			}
			else{
				throw new FormatException("Error due to : "+f);
			}
		}

		String getProperty(String k){
			try{
				return this.filtermap.get(k);
			}
			catch(Exception e){
				System.err.println("Key does not exist : "+k);
				return "";
			}
		}
	}

	public class FormatField{
		Map<String, String> formatmap=new LinkedHashMap<String, String>();
		public FormatField(){}
		public FormatField(String f) throws FormatException{
			if(f.startsWith("##FORMAT")){
				f=f.substring(10,f.length()-1);
				String[] sf=f.split(",");
				for(String f1:sf){
					String[] sf1=f1.split("=");
					if(sf1.length>2){ throw new FormatException("Error due to : "+f1);}
					try{
						this.formatmap.put(sf[0].trim(),sf[1].trim());
					}
					catch(Exception e){
						throw new FormatException("Error due to : "+f1);
					}
				}
			}
			else{
				throw new FormatException("Error due to : "+f);
			}
		}
		String getProperty(String k){
			try{
				return this.formatmap.get(k);
			}
			catch(Exception e){
				System.err.println("Key does not exist : "+k);
				return "";
			}
		}
	}

	public class AltField{
		Map<String, String> altmap=new LinkedHashMap<String, String>();
		public AltField(){}
		public AltField(String f) throws FormatException{
			if(f.startsWith("##ALT")){
				f=f.substring(7,f.length()-1);
				String[] sf=f.split(",");
				for(String f1:sf){
					String[] sf1=f1.split("=");
					if(sf1.length>2){ throw new FormatException("Error due to : "+f1);}
					try{
						this.altmap.put(sf[0].trim(),sf[1].trim());
					}
					catch(Exception e){
						throw new FormatException("Error due to : "+f1);
					}
				}
			}
			else{
				throw new FormatException("Error due to : "+f);
			}
		}
		String getProperty(String k){
			try{
				return this.altmap.get(k);
			}
			catch(Exception e){
				System.err.println("Key does not exist : "+k);
				return "";
			}
		}
	}

	public class AssemblyField{
		String assembly="";
		public AssemblyField(){}
		public AssemblyField(String f) throws FormatException{
			if(f.startsWith("##assembly")){
				assembly=f.substring(11).trim();
			}
			else{
				throw new FormatException("Error due to : "+f);
			}
		}
		String getAssembly(){
			return this.assembly;
		}
	}


	public class ReferenceField{
		String reference="";
		public ReferenceField(){}
		public ReferenceField(String f) throws FormatException{
			if(f.startsWith("##reference")){
				reference=f.substring(12).trim();
			}
			else{
				throw new FormatException("Error due to : "+f);
			}
		}
		String getReference(){
			return this.reference;
		}
	}

	public class ContigField{
		Map<String, String> contigmap=new LinkedHashMap<String, String>();
		public ContigField(){}
		public ContigField(String f) throws FormatException{
			if(f.startsWith("##contig")){
				f=f.substring(10,f.length()-1);
				String[] sf=f.split(",");
				for(String f1:sf){
					String[] sf1=f1.split("=");
					if(sf1.length>2){ throw new FormatException("Error due to : "+f1);}
					try{
						this.contigmap.put(sf[0].trim(),sf[1].trim());
					}
					catch(Exception e){
						throw new FormatException("Error due to : "+f1);
					}
				}
			}
			else{
				throw new FormatException("Error due to : "+f);
			}
		}
		String getProperty(String k){
			try{
				return this.contigmap.get(k);
			}
			catch(Exception e){
				System.err.println("Key does not exist : "+k);
				return "";
			}
		}
	}

	public class PedigreeField{
		Map<String, String> pedigreemap=new LinkedHashMap<String, String>();
		public PedigreeField(){}
		public PedigreeField(String f) throws FormatException{
			if(f.startsWith("##PEDIGREE")){
				f=f.substring(12,f.length()-1);
				String[] sf=f.split(",");
				for(String f1:sf){
					String[] sf1=f1.split("=");
					if(sf1.length>2){ throw new FormatException("Error due to : "+f1);}
					try{
						this.pedigreemap.put(sf[0].trim(),sf[1].trim());
					}
					catch(Exception e){
						throw new FormatException("Error due to : "+f1);
					}
				}
			}
			else{
				throw new FormatException("Error due to : "+f);
			}
		}
		String getProperty(String k){
			try{
				return this.pedigreemap.get(k);
			}
			catch(Exception e){
				System.err.println("Key does not exist : "+k);
				return "";
			}
		}
	}

	public class SampleField{
		Map<String, String> samplemap=new LinkedHashMap<String, String>();
		public SampleField(){}
		public SampleField(String f) throws FormatException{
			if(f.startsWith("##SAMPLE")){
				f=f.substring(10,f.length()-1);
				String[] sf=f.split(",");
				for(String f1:sf){
					String[] sf1=f1.split("=");
					if(sf1.length>2){ throw new FormatException("Error due to : "+f1);}
					try{
						this.samplemap.put(sf[0].trim(),sf[1].trim());
					}
					catch(Exception e){
						throw new FormatException("Error due to : "+f1);
					}
				}
			}
			else{
				throw new FormatException("Error due to : "+f);
			}
		}
		String getProperty(String k){
			try{
				return this.samplemap.get(k);
			}
			catch(Exception e){
				System.err.println("Key does not exist : "+k);
				return "";
			}
		}
	}

	public class MiscField{
		List<String> miscpair=new ArrayList<String>();
		//List containing key at 0 and value at 1
		public MiscField(){}
		public MiscField(String l) throws FormatException {
			if(l.startsWith("##")){
				int indeq=l.indexOf('=');
				String k=l.substring(2,indeq-1);
				String v=l.substring(indeq+1);
				miscpair.add(k);
				miscpair.add(v);
			}
			else{
				throw new FormatException("Error due to : ");
			}
		}
		String getProperty(String k){
			String prop="";
			if(miscpair.get(0)==k){
				prop=miscpair.get(1);
			}
			return prop;
		}
	}
	/*
	 * A class containing all the meta data class objects in a VCF file
	 */

	public class MetaData{
		List<String> fields=new ArrayList<String>();
		List<InfoField> infos=new ArrayList<InfoField>();
		List<FilterField> filters=new ArrayList<FilterField>();
		List<AltField> alts=new ArrayList<AltField>();
		List<ContigField> contigs=new ArrayList<ContigField>();
		List<FormatField> formats=new ArrayList<FormatField>();
		List<SampleField> samples=new ArrayList<SampleField>();
		List<PedigreeField> pedigrees=new ArrayList<PedigreeField>();
		Header header=new Header();
		AssemblyField assembly=new AssemblyField();
		ReferenceField reference=new ReferenceField();
		SourceField source=new SourceField();
		FileFormatField fileformat=new FileFormatField();

		List<MiscField> miscs=new ArrayList<MiscField>();

		int metafindexend=0;

		public MetaData(){
			//monofields
			fields.add("fileformat");	//0
			fields.add("source");
			fields.add("assembly");
			fields.add("reference");
			fields.add("contig");
			//multifields
			fields.add("INFO");		//5
			fields.add("FILTER");
			fields.add("FORMAT");
			fields.add("ALT");
			fields.add("PEDIGREE");
			fields.add("SAMPLE");		//10
		}
		//Mono-entities - entities present only once in the VCF file
		//fileformat, assembly, header, source, reference
		//Poly-entities - entities present multiple times in the VCF file
		//info, contig, alt, filter, format, sample
		//public MetaData(String fname) throws IOException, FormatException{
		public MetaData(VCFReader vc) throws IOException, FormatException{
			//monofields
			fields.add("fileformat");	//0
			fields.add("source");
			fields.add("assembly");
			fields.add("reference");
			//multifields
			fields.add("contig");
			fields.add("INFO");		//5
			fields.add("FILTER");
			fields.add("FORMAT");
			fields.add("ALT");
			fields.add("PEDIGREE");
			fields.add("SAMPLE");		//10
			String line;
			//BufferedReader br=new BufferedReader(new InputStreamReader(new FileInputStream(fname)));
			int i=0;
			while( (line=vc.readLine()) != null ){
				if(line.startsWith("#CHROM")){
					this.header=new Header(line);
					i++;
					break;
				}
				else if(line.startsWith("##")) //All the metadata conditions to be checked
				{
					i++;
					String field=line.substring(2,line.indexOf("=")-1);
					int fieldindex=fields.indexOf(field);
					switch(fieldindex){
						case -1:
							miscs.add(new MiscField(line));
							break;
						case 0:
							fileformat=new FileFormatField(line);
							break;
						case 1:
							source=new SourceField(line);
							break;
						case 2:
							assembly=new AssemblyField(line);
							break;
						case 3:
							reference=new ReferenceField(line);
							break;
						case 4:
							contigs.add(new ContigField(line));
							break;
						case 5:
							infos.add(new InfoField(line));
							break;
						case 6:
							filters.add(new FilterField(line));
							break;
						case 7:
							formats.add(new FormatField(line));
							break;
						case 8:
							alts.add(new AltField(line));
							break;
						case 9:
							pedigrees.add(new PedigreeField(line));
							break;
						case 10:
							samples.add(new SampleField(line));
							break;
					}
				}
			}
			metafindexend=i;
		}
		int variantStartLineNumber(){
			return this.metafindexend;
		}
		boolean hasInfoField(String fi){
			boolean hasinfo=false;
			for(int i=0; i < this.infos.size(); i++)
			{
				if(this.infos.get(i).getProperty("ID") == fi ){
					hasinfo=true;
				}
			}
			return hasinfo;
		}
	}

	public static class ANN{
		List<String> ann=new ArrayList<String>();
		public ANN(String s) throws FormatException{
			if(s.endsWith("|")){
				s+='.';
			}
			ann=Arrays.asList(s.split(Pattern.quote("|")));
		}
		public ANN(){}
		String getAllele(){
			return ann.get(0);
		}
		String getEffect(){
			return ann.get(1);
		}
		String getImpact(){
			return ann.get(2);
		}
		String getGene(){
			return ann.get(3);
		}
		String getGeneId(){
			return ann.get(4);
		}
		String getFeature(){
			return ann.get(5);
		}
		String getFeatureId(){
			return ann.get(6);
		}
		String getBiotype(){
			return ann.get(7);
		}
		String getRank(){
			return ann.get(8);
		}
		String getHGVSC(){
			return ann.get(9);
		}
		String getHGVSP(){
			return ann.get(10);
		}
		String getCDNAPos(){
			return ann.get(11);
		}
		String getCDNALen(){
			return ann.get(12);
		}
		String getCDSPos(){
			return ann.get(13);
		}
		String getCDSLen(){
			return ann.get(14);
		}
		String getAAPos(){
			return ann.get(15);
		}
		String getAALen(){
			return ann.get(16);
		}
		String getDistance(){
			return ann.get(17);
		}
		String getErrors(){
			return ann.get(18);
		}
	}

	public static class CSQ{
		List<String> csq=new ArrayList<String>();
		//single CSQ passed to the constructor
		public CSQ(String s) throws FormatException{
			if(s.endsWith("|")){
				s+=".";
			}
			csq=Arrays.asList(s.split(Pattern.quote("|")));
			if(csq.size()!=59){
				throw new FormatException("Error due to : "+s);
			}
		}
		public CSQ(){}//Start CSQ field extraction//
		String getAllele(){
			return this.csq.get(0);
		}
		String getConsequence(){
			return this.csq.get(1);
		}
		String getImpact(){
			return this.csq.get(2);
		}
		String getSymbol(){
			return this.csq.get(3);
		}
		String getGene(){
			return this.csq.get(4);
		}
		String getFeatureType(){
			return this.csq.get(5);
		}
		String getFeature(){
			return this.csq.get(6);
		}
		String getBiotype(){
			return this.csq.get(7);
		}
		String getExon(){
			return this.csq.get(8);
		}
		String getIntron(){
			return this.csq.get(9);
		}
		String getHGVSc(){
			return this.csq.get(10);
		}
		String getHGVSp(){
			return this.csq.get(11);
		}
		String getCDNAPosition(){
			return this.csq.get(12);
		}
		String getCDSPosition(){
			return this.csq.get(13);
		}
		String getProteinPosition(){
			return this.csq.get(14);
		}
		String getAminoAcids(){
			return this.csq.get(15);
		}
		String getCodons(){
			return this.csq.get(16);
		}
		String getExistingVariation(){
			return this.csq.get(17);
		}
		String getAlleleNum(){
			return this.csq.get(18);
		}
		String getDistance(){
			return this.csq.get(19);
		}
		String getStrand(){
			return this.csq.get(20);
		}
		String getVariantClass(){
			return this.csq.get(21);
		}
		String getMinimised(){
			return this.csq.get(22);
		}
		String getSymbolSource(){
			return this.csq.get(23);
		}
		String getHGNCId(){
			return this.csq.get(24);
		}
		String getCanonical(){
			return this.csq.get(25);
		}
		String getTSL(){
			return this.csq.get(26);
		}
		String getCCDS(){
			return this.csq.get(27);
		}
		String getENSP(){
			return this.csq.get(28);
		}
		String getSwissprot(){
			return this.csq.get(29);
		}
		String getTREMBL(){
			return this.csq.get(30);
		}
		String getUniparc(){
			return this.csq.get(31);
		}
		String getSIFT(){
			return this.csq.get(32);
		}
		String getPolyPhen(){
			return this.csq.get(33);
		}
		String getDomains(){
			return this.csq.get(34);
		}
		String getHGSVOffset(){
			return this.csq.get(35);
		}
		String getGMaf(){
			return this.csq.get(36);
		}
		String getAfrMaf(){
			return this.csq.get(37);
		}
		String getAmrMaf(){
			return this.csq.get(38);
		}
		String getAsnMaf(){
			return this.csq.get(39);
		}
		String getEasMaf(){
			return this.csq.get(40);
		}
		String getEurMaf(){
			return this.csq.get(41);
		}
		String getSasMaf(){
			return this.csq.get(42);
		}
		String getAaMaf(){
			return this.csq.get(43);
		}
		String getEaMaf(){
			return this.csq.get(44);
		}
		String getClinSig(){
			return this.csq.get(45);
		}
		String getSomatic(){
			return this.csq.get(46);
		}
		String getPheno(){
			return this.csq.get(47);
		}
		String getPubmed(){
			return this.csq.get(48);
		}
		String getMotifName(){
			return this.csq.get(49);
		}
		String getMotifPos(){
			return this.csq.get(50);
		}
		String getHighInfPos(){
			return this.csq.get(51);
		}
		String getMotifScoreChange(){
			return this.csq.get(52);
		}
		String getLofInfo(){
			return this.csq.get(53);
		}
		String getLofFlags(){
			return this.csq.get(54);
		}
		String getLofFilter(){
			return this.csq.get(55);
		}
		String getLof(){
			return this.csq.get(56);
		}
		String getContext(){
			return this.csq.get(57);
		}
		String getAncestral(){
			return this.csq.get(58);
		}
	}


	/* Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|ALLELE_NUM|DISTANCE|STRAND|VARIANT_CLASS|MINIMISED|SYMBOL_SOURCE|HGNC_ID|CANONICAL|TSL|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|SIFT|PolyPhen|DOMAINS|HGVS_OFFSET|GMAF|AFR_MAF|AMR_MAF|ASN_MAF|EAS_MAF|EUR_MAF|SAS_MAF|AA_MAF|EA_MAF|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|LoF_info|LoF_flags|LoF_filter|LoF|context|ancestral * */
	public class InfoCSQ{
		//may contain multiple CSQ, comma separated - functional concequences
		List<CSQ> csqs=new ArrayList<CSQ>();

		public InfoCSQ(){}
		public InfoCSQ(String ic) throws FormatException {
			for(String t:ic.split(",")){
				this.csqs.add(new CSQ(t));
			}
		}
		int getNumberOfCSQs(){
			return this.csqs.size();
		}

		CSQ getCSQ(int i){
			return this.csqs.get(i);
		}
	}
	public class InfoANN{
		List<ANN> anns=new ArrayList<ANN>();
		public InfoANN(){}
		public InfoANN(String an) throws FormatException{
			for(String t:an.split(",")){
				this.anns.add(new ANN(t));
			}
		}
		int getNumberOfANNs(){
			return this.anns.size();
		}
		ANN getANN(int i){
			return this.anns.get(i);
		}
	}
	/*
	 * Variant entry class
	 */

	public static class Variant{
		//String varentry="";
		String chrom="", ref="", filter="", qual="";
		int pos;
		List<String> id = new ArrayList<String>(), alt = new ArrayList<String>(), additional_fields=new ArrayList<String>();
		Map<String,String> info=new LinkedHashMap<String, String>();
		public Variant(){}
		public Variant(String v) throws FormatException {
			//varentry=v;
			List<String> sv=new ArrayList<String>();
			for(String tmp:v.trim().split("\t"))
			{
				sv.add(tmp);
			}
			//checks on sv length
			if(sv.size()<8){
				throw new FormatException("Error due to : "+v);
			}
			else{
				this.chrom=sv.get(0);
				this.pos=Integer.parseInt(sv.get(1));
				if(sv.get(2).contains(";")){
					for(String tmp:sv.get(2).split(";")){
							this.id.add(tmp);
					}
				}
				else if (sv.get(2).contains(",")) {
					for(String tmp:sv.get(2).split(",")){
							this.id.add(tmp);
					}
				}
				else{
						this.id.add(sv.get(2));
				}
				this.ref=sv.get(3);
				if(sv.get(4).contains(",")){ //ALT column
					for(String tmp:sv.get(4).split(",")){
						this.alt.add(tmp);
					}
				}
				else{
					this.alt.add(sv.get(4));
				}
				this.qual=sv.get(5);
				this.filter=sv.get(6);
				if(sv.get(7).contains(";")){
					for(String tmp:sv.get(7).split(";")){
						if(tmp.contains("=")){
							String[] arr=tmp.split("=");
							this.info.put(arr[0],arr[1]);
						}
						else{
							this.info.put(tmp,"");//For handling flags, which may not have a key->value pair
							//continue;
							//throw new FormatException("Error due to : "+v);
						}
					}
				}
				if (sv.size()>8) {
					for (int i=8; i < sv.size(); i++) {
						this.additional_fields.add(sv.get(i));
					}
				}
			}
		}

		public Variant(ArrayList<String> sv) throws FormatException {
			//checks on sv length
			if(sv.size()<8){
				throw new FormatException("Error: Length of list is less than 8 (expected size of VCF entry)!!");
			}
			else{
				this.chrom=sv.get(0);
				this.pos=Integer.parseInt(sv.get(1));
				if(sv.get(2).contains(";")){
					for(String tmp:sv.get(2).split(";")){
							this.id.add(tmp);
					}
				}
				else if (sv.get(2).contains(",")) {
					for(String tmp:sv.get(2).split(",")){
							this.id.add(tmp);
					}
				}
				else{
						this.id.add(sv.get(2));
				}
				this.ref=sv.get(3);
				if(sv.get(4).contains(",")){ //ALT column
					for(String tmp:sv.get(4).split(",")){
						this.alt.add(tmp);
					}
				}
				else{
					this.alt.add(sv.get(4));
				}
				this.qual=sv.get(5);
				this.filter=sv.get(6);
				if(sv.get(7).contains(";")){
					for(String tmp:sv.get(7).split(";")){
						if(tmp.contains("=")){
							String[] arr=tmp.split("=");
							this.info.put(arr[0],arr[1]);
						}
						else{
							this.info.put(tmp,"");//For handling flags, which may not have a key->value pair
							//continue;
							//throw new FormatException("Error due to : "+v);
						}
					}
				}
				if (sv.size()>8) {
					for (int i=8; i < sv.size(); i++) {
						this.additional_fields.add(sv.get(i));
					}
				}
			}
		}

		String getChromosome(){
			return this.chrom;
		}
		List<String> getVariantKeys(){
			List<String> vkeys=new ArrayList<String>();
			for(String aa:this.getAlternateAlleles()){
				vkeys.add(this.getChromosome()+"#"+Integer.toString(this.getPosition())+"#"+this.getReferenceAllele()+"#"+aa); //variant key
			}
			return vkeys;
		}
		boolean isSNP(){
			boolean issnp=false;
			for (String aa : this.alt) {
					if (aa.length()==1) {
						issnp=true;
					}
			}
			return this.ref.length()==1 && issnp;
		}
		private boolean isBasePurine(String base){
				boolean pur=false;
				if (base.equals("A") || base.equals("G")) {
					pur=true;
				}
				return pur;
		}
		private boolean isBasePyrimidine(String base){
			boolean pyr=false;
			if (base.equals("C") || base.equals("T")) {
				pyr=true;
			}
			return pyr;
		}
		boolean isTransition(){
			boolean istra=false;
			if (this.isSNP()) {
				for (String aa: this.alt ) {
					if (this.isBasePurine(this.ref) && this.isBasePurine(aa)) {
						istra=true;
						break;
					}
					else if (this.isBasePyrimidine(this.ref) && this.isBasePyrimidine(aa)) {
						istra=true;
						break;
					}
				}
			}
			return istra;
		}
		boolean isTransversion(){
			boolean istra=false;
			if (this.isSNP()) {
				for (String aa: this.alt ) {
					if (this.isBasePurine(this.ref) && this.isBasePyrimidine(aa)) {
						istra=true;
						break;
					}
					else if (this.isBasePyrimidine(this.ref) && this.isBasePurine(aa)) {
						istra=true;
						break;
					}
				}
			}
			return istra;
		}
		boolean hasMultipleAlleles(){
			return this.alt.size()>1;
		}
		int getNumberOfAlternateAlleles(){
			return this.alt.size();
		}
		List<String> getAlternateAlleles(){
			return this.alt;
		}
		void setAlternateAlleles(List<String> aa){
			this.alt=aa;
		}
		GenomicInterval toGenomicInterval() throws GenomicInterval.GenomicIntervalException {
			return new GenomicInterval(this.getChromosome(), this.getPosition(),this.getPosition()+this.ref.length()-1);
		}

		boolean equals(Variant v){
			boolean areequal=false;
			if(this.chrom.equals(v.getChromosome()) && this.pos==v.getPosition() && this.ref.equals(v.getReferenceAllele())){
				for(String al : this.id){
					for(String val : v.getIds()){
						if(al.equals(val)){
							areequal=true;
							break;
						}
					}
				}
			}
			return areequal;
		}

		String alternateAlleleToString(){
			return String.join(",",this.alt);
		}
		String getAlternateAllele(int i) throws IndexOutOfBoundsException {
			return this.alt.get(i);
		}
		void setAlternateAllele(int i, String e) throws IndexOutOfBoundsException {
			this.alt.set(i, e);
		}
		String getReferenceAllele(){
			return this.ref;
		}
		void setReferenceAllele(String r){
			this.ref=r;
		}
		int getPosition(){
			return this.pos;
		}
		void setPosition(int i){
			this.pos=i;
		}
		String getFilter(){
			return this.filter;
		}
		String getQuality(){
			return this.qual;
		}
		String getInfo(String k){
			try{
				return this.info.get(k);
			}
			catch(Exception e){
				System.err.println("Key not found - "+k);
				return "";
			}
		}
		String infoToString(){
			String sinfo="", resinfo="";
			for(Map.Entry<String,String> entry : this.info.entrySet())
			{
				if(entry.getValue().equals(""))
				{
					sinfo=sinfo+entry.getKey()+";";
				}
				else{
					sinfo=sinfo+entry.getKey()+"="+entry.getValue()+";";
				}
			}
			if (sinfo.length()>0) {
				resinfo=sinfo.substring(0, sinfo.length()-1);
			}
			return resinfo;
		}
		boolean hasInfo(String k){
			return info.containsKey(k);
		}
		boolean isInfoFlag(String k){
			boolean isflag=false;
			if(this.hasInfo(k)){
				if(this.getInfo(k)==""){
					isflag=true;
				}
			}
			else{
				System.err.println("Key not found - "+k);
			}
			return isflag;
		}
		String getId(int i){
			try{
				return this.id.get(i);
			}
			catch(Exception e){
				System.err.println("Out of index!!");
				return "";
			}
		}
		void addId(String ii){
			this.id.add(ii);
		}
		String idToString(){
			String recid= String.join(";",this.id);
			if(recid.startsWith(".;")){
				recid=recid.substring(2,recid.length());
			}
			return recid;
		}
		List<String> getIds(){
			return this.id;
		}
		String getIdsAsString(){
			String ids="";
			for (int g=0; g < this.id.size() ; g++ ) {
				if(ids.length()>1){
					ids=ids+","+this.id.get(g);
				}
				else{
					ids=this.id.get(g);
				}
			}
			return ids;
		}
		int getNumberOfIds(){
			return this.id.size();
		}
		public String toString(){
			String rec=this.getChromosome()+"\t"+Long.toString(this.pos)+"\t"+this.idToString()+"\t"+this.ref+"\t"+this.alternateAlleleToString()+"\t"+this.qual+"\t"+this.filter+"\t"+this.infoToString();
			if (additional_fields.size()>0) {
				for (int i=0; i<additional_fields.size() ; i++ ) {
					rec=rec+"\t"+additional_fields.get(i);
				}
			}
			return rec;
		}
		ArrayList<VCFReader.Variant> getAlternateAlleleAsVariants() throws FormatException {
			ArrayList<VCFReader.Variant> altvars=new ArrayList<VCFReader.Variant>();
			for(int a = 0 ; a < this.alt.size(); a++){
				String rec=this.getChromosome()+"\t"+Long.toString(this.pos)+"\t"+this.idToString()+"\t"+this.ref+"\t"+this.getAlternateAllele(a)+"\t"+this.qual+"\t"+this.filter+"\t"+this.infoToString();
				if (additional_fields.size()>0) {
					for (int i=0; i<additional_fields.size() ; i++ ) {
						rec=rec+"\t"+additional_fields.get(i);
					}
				}
				altvars.add(new VCFReader.Variant(rec));
			}
			return altvars;
		}
		List<CSQ> getCSQs() throws FormatException{
			List<CSQ> css=new ArrayList<CSQ>();
			if(this.hasInfo("CSQ")){
				String inf=this.getInfo("CSQ");
				for(String t:inf.split(",")){
					css.add(new CSQ(t));
				}
			}
			return css;
		}
		CSQ getCSQ(int i) throws FormatException {
			return this.getCSQs().get(i);
		}
		int getNumberOfCSQs() throws FormatException{
			return this.getCSQs().size();
		}
		List<ANN> getANNs() throws FormatException{
			List<ANN> ann=new ArrayList<ANN>();
			if(this.hasInfo("ANN")){
				String inf=this.getInfo("ANN");
				for(String t:inf.split(",")){
					ann.add(new ANN(t));
				}
			}
			else{
				System.out.println("Record does not contain ANN field!!");
			}
			return ann;
		}
		ANN getANN(int i) throws FormatException {
			return this.getANNs().get(i);
		}
		int getNumberOfANNs() throws FormatException{
			return this.getANNs().size();
		}
	}
	//Test for alternate allele variant separation - Sep 4, 2019 - pass
	/*
	public static void main(String[] args){
		try{
			System.out.println(args[0]);
			VCFReader v=new VCFReader(args[0]);
			Variant var=v.getNextVariant();
			List<Variant> avars=var.getAlternateAlleleAsVariants();
			for(int iv=0; iv < avars.size(); iv++){
				System.out.println(avars.get(iv).toString());
			}
		}
		catch(Exception e){
			System.out.println(args[0]);
			e.printStackTrace();
		}
	}

	//Test for working of VCFReader class
	public static void main(String[] args){
		try{
			TabixReader tr=new TabixReader(args[1]);
			String s;
			VCFReader v=new VCFReader(args[0]);
			while(true)
			{
				try{
					Variant var=v.getNextVariant();
					for(int i=0; i < var.getNumberOfCSQs(); i++){
					//	System.out.println(var.getChromosome()+" "+var.getPosition()+" "+var.getCSQ(i).getConsequence());
						TabixReader.Iterator iter=tr.query(var.chrom+":"+var.getPosition()+"-"+var.getPosition());
						while(iter!= null && (s=iter.next())!=null){
							System.out.println(s);
						}
					}
				}
				catch(NullPointerException ne){
					break;
				}
				catch(Exception e)
				{
					e.printStackTrace();
					break;
				}
			}
		}
		catch(Exception e){
			e.printStackTrace();
		}
	}
	*/
}
