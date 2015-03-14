/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package discord;

import JDDiscord.JDDiscord;
import edu.hawaii.jmotif.sax.SAXFactory;
import edu.hawaii.jmotif.sax.datastructures.DiscordRecords;
import weka.core.Instances;
import weka.core.converters.ConverterUtils.DataSource;

/**
 *
 * @author ian
 */
public class Discord {

    private static final int m =-1;

    // ERP
 //   private static final int SLIDING_WINDOW_SIZE = 128;
 //   private static final int ALPHABET_SIZE = 3;
 //   private static final int NUM_TO_REPORT = 1;
 //   private static final int J = 6;
 //   private static final String file = "erp\\erp.arff";
 //   private static final String DATA_VALUE_ATTRIBUTE = "value0";
    
    // tickwise
//    private static final int SLIDING_WINDOW_SIZE = 128;
//    private static final int ALPHABET_SIZE = 3;
//    private static final int NUM_TO_REPORT = 1;
//    private static final int J = 4;
//    private static final String file = "tickwise\\tickwise.arff";
//    private static final String DATA_VALUE_ATTRIBUTE = "value0";
    
    // koski
//    private static final int SLIDING_WINDOW_SIZE = 128;
//    private static final int ALPHABET_SIZE = 3;
//    private static final int NUM_TO_REPORT = 1;
//    private static final int J = 6;
//   private static final String file = "koski\\koski.arff";
//    private static final String DATA_VALUE_ATTRIBUTE = "value0";
    
    // qtdbsel102
//    private static final int SLIDING_WINDOW_SIZE = 128;
//    private static final int ALPHABET_SIZE = 3;
//    private static final int NUM_TO_REPORT = 1;
//    private static final int J = 6;
//    private static final String file = "qtdbsel102\\qtdbsel102.arff";
//    private static final String DATA_VALUE_ATTRIBUTE = "value0";
    
    //har
//    private static final int SLIDING_WINDOW_SIZE = 520;
//    private static final int ALPHABET_SIZE = 5;
//    private static final int NUM_TO_REPORT = 10;
//    private static final int J = 2;
//    private static final String file = "har\\activity01_30000-39999.arff";
//    private static final String DATA_VALUE_ATTRIBUTE = "value0";

    //hbgen
//    private static final int SLIDING_WINDOW_SIZE = 16;
//    private static final int ALPHABET_SIZE = 5;
//    private static final int NUM_TO_REPORT = 10;
//    private static final int J = 4;
//    private static final String file = "hbgen\\hb.arff";
//    private static final String DATA_VALUE_ATTRIBUTE = "diff";
    
    //exnoise
//    private static final int SLIDING_WINDOW_SIZE = 40;
//    private static final int ALPHABET_SIZE = 5;
//    private static final int NUM_TO_REPORT = 5;
//    private static final int J = 5;
//    private static final String file = "exnoise\\exnoise.arff";
//    private static final String DATA_VALUE_ATTRIBUTE = "value0";
    
    //ecg
//    private static final int SLIDING_WINDOW_SIZE = 180;
//    private static final int ALPHABET_SIZE = 5;
//    private static final int NUM_TO_REPORT = 10;
//    private static final int J = 5;
//    private static final String file = "ecg\\sel102s300_400.arff";
//    private static final String DATA_VALUE_ATTRIBUTE = "value0";

        //ecg102
    private static final int SLIDING_WINDOW_SIZE = 360;
    private static final int ALPHABET_SIZE = 5;
    private static final int NUM_TO_REPORT = 20;
    private static final int J = 1;
    private static final String file = "..\\parallel_discord_rank\\ecg102.arff";
    private static final String DATA_VALUE_ATTRIBUTE = "value0";
    
    /**
     * Executable method.
     *
     * @param args None used.
     * @throws Exception if error occurs.
     */
    public static void main(String[] args) throws Exception {

        // get the data first
        Instances tsData = readTSData(file);
        if (m != -1) {
            while (tsData.numInstances() > m) {
                tsData.delete(tsData.numInstances()-1);
            }
        }

        // now build the SAX data structure using sliding window of size 40 and alphabet of 3
        DiscordRecords jddr = SAXFactory.instances2Discords(tsData, DATA_VALUE_ATTRIBUTE, SLIDING_WINDOW_SIZE, ALPHABET_SIZE, NUM_TO_REPORT);
//        DiscordRecords jddr = JDDiscord.instanceJDDiscords(tsData, DATA_VALUE_ATTRIBUTE,
//                SLIDING_WINDOW_SIZE, ALPHABET_SIZE, NUM_TO_REPORT, J);

//        DiscordRecords jddr = JDDiscord.instanceLOFDiscords(tsData, DATA_VALUE_ATTRIBUTE,
//                SLIDING_WINDOW_SIZE, NUM_TO_REPORT, J);
//        Attribute dataAttribute = tsData.attribute(DATA_VALUE_ATTRIBUTE);
//        DiscordRecords jddr = JDDiscord.findBFDiscords(JDDiscord.toRealSeries(tsData, dataAttribute), SLIDING_WINDOW_SIZE,5);
        
        // printout the discords occurrences
//        System.out.print("\nFive best discords:\n\n" + jddr.toString());
//
//        System.out.print("\n\n\ndiscords=t(rev(c(" + jddr.toCoordinates()
//                + ")))\nwords=t(rev(c(" + jddr.toPayloads() + ")))\ndistances=t(rev(c(" + jddr.toDistances()
//                + ")))\n\nDONE");

    }

    /**
     * Read the timeseries data into WEKA format.
     *
     * @return Timeseries.
     * @throws Exception If error occurs.
     */
    private static Instances readTSData(String filename) throws Exception {
        Instances data = DataSource.read(filename);
        return data;
    }

}
