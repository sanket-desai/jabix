//Extracts missense using CSQ from ExAc
package jabix;
//arguments: vcffile, outputfile
import java.io.PrintWriter;
import java.time.Instant;
import java.time.Duration;
public class VCFExtractMissense{
	public static void main(String[] args){
		long variantsprocessed=0;
		Instant start=Instant.now();
		try{
			PrintWriter writer=new PrintWriter(args[1]);
			writer.println("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
			VCFReader v=new VCFReader(args[0]);
			while(true){
				try{
					VCFReader.Variant vr=v.getNextVariant();
					variantsprocessed++;
					for(int i=0; i<vr.getNumberOfCSQs(); i++){
						VCFReader.CSQ vrcs=vr.getCSQ(i);
						if(vrcs.getConsequence().startsWith("missense_variant") && vrcs.getSwissprot().length()>1){
							//System.out.println(vr.getChromosome()+"\t"+vr.getPosition()+"\t"+vr.getId(0)+"\t"+vr.getReferenceAllele()+"\t"+vr.getAlternateAllele(0)+"\t"+vr.getQuality()+"\t"+vr.getFilter()+"\t"+ "AF="+vr.getInfo("AF")+";Symbol="+vrcs.getSymbol()+";ProteinPosition="+vrcs.getProteinPosition()+";AAChange="+vrcs.getAminoAcids()+";Swissprot="+vrcs.getSwissprot()+";HGVSp="+vrcs.getHGVSp());
							writer.println(vr.getChromosome()+"\t"+vr.getPosition()+"\t"+vr.getId(0)+"\t"+vr.getReferenceAllele()+"\t"+vr.getAlternateAllele(0)+"\t"+vr.getQuality()+"\t"+vr.getFilter()+"\t"+ "AF="+vr.getInfo("AF")+";Symbol="+vrcs.getSymbol()+";ProteinPosition="+vrcs.getProteinPosition()+";AAChange="+vrcs.getAminoAcids()+";Swissprot="+vrcs.getSwissprot()+";HGVSp="+vrcs.getHGVSp());
						}
					}
				}
				catch(NullPointerException ne){
					break;
				}
				catch(Exception e){
					e.printStackTrace();
					break;
				}
			}
			writer.close();
		}
		catch(Exception e){
			e.printStackTrace();
		}
		Instant end=Instant.now();
		Duration timespent=Duration.between(start,end);
		System.out.println("Variants processed: "+variantsprocessed);
		System.out.println("Time taken: "+timespent.toHours()+ " hours!!");
		System.out.println("Time taken: "+timespent.toMinutes()+ " minutes!!");
		System.out.println("Time taken: "+timespent.toMillis()+ " milliseconds!!");
	}
}
