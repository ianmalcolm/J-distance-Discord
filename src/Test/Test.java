/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package Test;

import JDDiscord.DataHandler;
import JDDiscord.Distance;
import JDDiscord.JDDiscord;
import SAXFactory.DiscordRecords;
import SAXFactory.SAXFactory;
import SAXFactory.TSUtils;
import java.util.Date;
import java.util.logging.ConsoleHandler;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import weka.core.Attribute;
import weka.core.Instances;
import weka.core.converters.ConverterUtils.DataSource;

/**
 *
 * @author ian
 */
public class Test {

//    private static final Logger logger = Logger.getLogger(Test.class.getName());
    static {
//        logger.setLevel(Level.ALL);
        ConsoleHandler ch = new ConsoleHandler();
        ch.setLevel(Level.ALL);
        Logger.getLogger("").addHandler(ch);
    }

    // ERP
    //   private static final int SLIDING_WINDOW_SIZE = 128;
    //   private static final int ALPHABET_SIZE = 3;
    //   private static final int NUM_TO_REPORT = 1;
    //   private static final int J = 6;
    //   private static final String file = "erp/erp.arff";
    //   private static final String DATA_VALUE_ATTRIBUTE = "value0";
    // tickwise
//    private static final int SLIDING_WINDOW_SIZE = 128;
//    private static final int ALPHABET_SIZE = 3;
//    private static final int NUM_TO_REPORT = 1;
//    private static final int J = 4;
//    private static final String file = "tickwise/tickwise.arff";
//    private static final String DATA_VALUE_ATTRIBUTE = "value0";
    // koski
//    private static final int SLIDING_WINDOW_SIZE = 128;
//    private static final int ALPHABET_SIZE = 3;
//    private static final int NUM_TO_REPORT = 1;
//    private static final int J = 6;
//   private static final String file = "koski/koski.arff";
//    private static final String DATA_VALUE_ATTRIBUTE = "value0";
    // qtdbsel102
//    private static final int SLIDING_WINDOW_SIZE = 128;
//    private static final int ALPHABET_SIZE = 3;
//    private static final int NUM_TO_REPORT = 1;
//    private static final int J = 6;
//    private static final String file = "qtdbsel102/qtdbsel102.arff";
//    private static final String DATA_VALUE_ATTRIBUTE = "value0";
    //har
//    private static final int SLIDING_WINDOW_SIZE = 520;
//    private static final int ALPHABET_SIZE = 5;
//    private static final int NUM_TO_REPORT = 10;
//    private static final int J = 2;
//    private static final String file = "har/activity01_30000-39999.arff";
//    private static final String DATA_VALUE_ATTRIBUTE = "value0";
    //hbgen
//    private static final int SLIDING_WINDOW_SIZE = 16;
//    private static final int ALPHABET_SIZE = 5;
//    private static final int NUM_TO_REPORT = 10;
//    private static final int J = 4;
//    private static final String file = "hbgen/hb.arff";
//    private static final String DATA_VALUE_ATTRIBUTE = "diff";
    //exnoise
//    private static final int SLIDING_WINDOW_SIZE = 40;
//    private static final int ALPHABET_SIZE = 5;
//    private static final int NUM_TO_REPORT = 5;
//    private static final int J = 5;
//    private static final String file = "exnoise/exnoise.arff";
//    private static final String DATA_VALUE_ATTRIBUTE = "value0";
    //ecg
//    private static final int SLIDING_WINDOW_SIZE = 180;
//    private static final int ALPHABET_SIZE = 5;
//    private static final int NUM_TO_REPORT = 10;
//    private static final int J = 5;
//    private static final String file = "ecg/sel102s300_400.arff";
//    private static final String DATA_VALUE_ATTRIBUTE = "value0";
    //ecg102
//    private static final int SLIDING_WINDOW_SIZE = 360;
//    private static final int ALPHABET_SIZE = 5;
//    private static final int NUM_TO_REPORT = 20;
//    private static final int J = 1;
//    private static final String file = "../parallel_discord_rank/ecg102.arff";
//    private static final String DATA_VALUE_ATTRIBUTE = "value0";
    /**
     * Executable method.
     *
     * @param args None used.
     * @throws Exception if error occurs.
     */
    public static void main(String[] args) throws Exception {

        Date totalstart = new Date();

        int SLIDING_WINDOW_SIZE = 360;
        int ALPHABET_SIZE = 5;
        String DATA_VALUE_ATTRIBUTE = "value0";
        String FILE = "../datasets/ecg/ecg102.arff";
        int DIMENSION = 5;
        int LENGTH = 8000;
        int REPORT_NUM = 5;
        boolean ENHANCED = true;
        int J = 3;
        Level level = Level.INFO;

        if (args.length > 0) {
            Options options = new Options();
            options.addOption("len", true, "Set the length of dataset, -1 for using the whole dataset");
            options.addOption("rep", true, "The number of reported discords");
            options.addOption("j", true, "the number of nearest neighbor get involed in calculation");
            options.addOption("fil", true, "The file name of the dataset");
            options.addOption("att", true, "The abbribute of instances, see the introduction of arff in WEKA for details");
            options.addOption("alp", true, "The size of alphabets");
            options.addOption("win", true, "The size of sliding window");
            options.addOption("dim", true, "The dimension of SAXTrie, currently it is internally set the same as the size of alphabets");
            options.addOption("log", true, "The log level");
            options.addOption("enh", false, "enhanced early abandon technique");
            options.addOption("h", false, "Print help message");

            CommandLineParser parser = new BasicParser();
            CommandLine cmd = parser.parse(options, args);

            if (cmd.hasOption("h")) {
                HelpFormatter formatter = new HelpFormatter();
                formatter.printHelp("jdd", options);
                return;
            }
            if (cmd.hasOption("enh")) {
                ENHANCED = true;
            }
            if (cmd.hasOption("alp")) {
                ALPHABET_SIZE = Integer.parseInt(cmd.getOptionValue("alp"));
            }
            if (cmd.hasOption("win")) {
                SLIDING_WINDOW_SIZE = Integer.parseInt(cmd.getOptionValue("win"));
            }
            if (cmd.hasOption("len")) {
                LENGTH = Integer.parseInt(cmd.getOptionValue("len"));
            }
            if (cmd.hasOption("rep")) {
                REPORT_NUM = Integer.parseInt(cmd.getOptionValue("rep"));
            }
            if (cmd.hasOption("j")) {
                J = Integer.parseInt(cmd.getOptionValue("j"));
            }

            if (cmd.hasOption("dim")) {
                DIMENSION = Integer.parseInt(cmd.getOptionValue("dim"));
            }
            if (cmd.hasOption("fil")) {
                FILE = cmd.getOptionValue("fil");
            }
            if (cmd.hasOption("att")) {
                DATA_VALUE_ATTRIBUTE = cmd.getOptionValue("att");
            }

            if (cmd.hasOption("log")) {
                if (cmd.getOptionValue("log").equalsIgnoreCase("OFF")) {
                    level = Level.OFF;
                } else if (cmd.getOptionValue("log").equalsIgnoreCase("SEVERE")) {
                    level = Level.SEVERE;
                } else if (cmd.getOptionValue("log").equalsIgnoreCase("WARNING")) {
                    level = Level.WARNING;
                } else if (cmd.getOptionValue("log").equalsIgnoreCase("INFO")) {
                    level = Level.INFO;
                } else if (cmd.getOptionValue("log").equalsIgnoreCase("CONFIG")) {
                    level = Level.CONFIG;
                } else if (cmd.getOptionValue("log").equalsIgnoreCase("FINE")) {
                    level = Level.FINE;
                } else if (cmd.getOptionValue("log").equalsIgnoreCase("FINER")) {
                    level = Level.FINER;
                } else if (cmd.getOptionValue("log").equalsIgnoreCase("FINEST")) {
                    level = Level.FINEST;
                }
            }
        }

        System.out.println("-fil " + FILE
                + " -att " + DATA_VALUE_ATTRIBUTE
                + " -win " + SLIDING_WINDOW_SIZE
                + " -len " + LENGTH
                + " -rep " + REPORT_NUM
                + " -alp " + ALPHABET_SIZE
                + " -dim " + DIMENSION
                + " -j " + J
                + " -log " + level.toString()
                + (ENHANCED ? " -enh" : ""));

        JDDiscord.setLoggerLevel(level);

        // get the data first
        Instances tsData = readTSData(FILE);
        if (LENGTH != -1) {
            while (tsData.numInstances() > LENGTH) {
                tsData.delete(tsData.numInstances() - 1);
            }
        }

        Attribute dataAttribute = tsData.attribute(DATA_VALUE_ATTRIBUTE);
        double[] series = SAXFactory.toRealSeries(tsData, dataAttribute);
        DataInMemory dh = new DataInMemory(series, SLIDING_WINDOW_SIZE);
        ED ed = new ED();
        JDDiscord jdd = new JDDiscord(series, SLIDING_WINDOW_SIZE, ALPHABET_SIZE, DIMENSION, dh, ed);
        jdd.setEnhanced(ENHANCED);

        // now build the SAX data structure using sliding window of size 40 and alphabet of 3
        Date jddstart = new Date();
        DiscordRecords jddr = jdd.findJDDiscords(REPORT_NUM, J);
        Date jddend = new Date();

        Date totalend = new Date();

        // printout the discords occurrences
        System.out.print("\nFive best discords:\n\n" + jddr.toString());
        System.out.println();

        System.out.println("Discord discovery elapsed time: " + SAXFactory.timeToString(jddstart.getTime(), jddend.getTime()));
        System.out.println("Total elapsed time: " + SAXFactory.timeToString(totalstart.getTime(), totalend.getTime()));
        System.out.println("Total count of the calls to the distance function: " + jdd.totalcnt);

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

class DataInMemory extends DataHandler {

    private double[] series;
    private int windowSize;
    private double mean;
    private double std;

    public DataInMemory(double[] _series, int _windowSize) {
        series = _series;
        windowSize = _windowSize;
        mean = TSUtils.mean(series);
        std = TSUtils.stDev(series);
    }

    @Override
    public long size() {
        return series.length - windowSize + 1;
    }

    @Override
    public double[] get(long i) {
        double[] subSeries = TSUtils.getSubSeries(series, (int) i, (int) i + windowSize);
        return TSUtils.zNormalize(subSeries, mean, std);
    }

}

class ED extends Distance {

    private long cnt = 0;

    /**
     * Calculates the square of the Euclidean distance between two 1D points
     * represented by real values.
     *
     * @param p1 The first point.
     * @param p2 The second point.
     * @return The Square of Euclidean distance.
     */
    public static double distance2(double p1, double p2) {
        double temp = p1 - p2;
        return temp * temp;
    }

    /**
     * Calculates the square of the Euclidean distance between two
     * multidimensional points represented by the real vectors.
     *
     * @param point1 The first point.
     * @param point2 The second point.
     * @return The Euclidean distance.
     * @throws TSException In the case of error.
     */
    public static double distance2(double[] point1, double[] point2) {
        assert point1.length == point2.length : "Exception in Euclidean distance: array lengths are not equal";
        Double sum = 0D;
        for (int i = 0; i < point1.length; i++) {
            double temp = point2[i] - point1[i];
            sum = sum + temp * temp;
        }
        return sum;
    }

    /**
     * Calculates the square of the Euclidean distance between two
     * multidimensional points represented by integer vectors.
     *
     * @param point1 The first point.
     * @param point2 The second point.
     * @return The Euclidean distance.
     * @throws TSException In the case of error.
     */
    public static double distance2(int[] point1, int[] point2) {
        assert point1.length == point2.length : "Exception in Euclidean distance: array lengths are not equal";
        Double sum = 0D;
        for (int i = 0; i < point1.length; i++) {
            double temp = Integer.valueOf(point2[i]).doubleValue() - Integer.valueOf(point1[i]).doubleValue();
            sum = sum + temp * temp;
        }
        return sum;
    }

    /**
     * Calculates the Euclidean distance between two points.
     *
     * @param p1 The first point.
     * @param p2 The second point.
     * @return The Euclidean distance.
     */
    public static double distance(double p1, double p2) {
        double temp = (p1 - p2);
        double d = temp * temp;
        return Math.sqrt(d);
    }

    /**
     * Calculates the Euclidean distance between two points.
     *
     * @param point1 The first point.
     * @param point2 The second point.
     * @return The Euclidean distance.
     * @throws TSException In the case of error.
     */
    @Override
    public double distance(double[] point1, double[] point2) {
        cnt++;
        return Math.sqrt(distance2(point1, point2));
//        return (distance2(point1, point2));
    }

    /**
     * Calculates the Euclidean distance between two points.
     *
     * @param point1 The first point.
     * @param point2 The second point.
     * @return The Euclidean distance.
     * @throws TSException In the case of error.
     */
    public static double distance(int[] point1, int[] point2) {
        return Math.sqrt(distance2(point1, point2));
    }

    /**
     * Calculates euclidean distance between two one-dimensional time-series of
     * equal length.
     *
     * @param series1 The first series.
     * @param series2 The second series.
     * @return The eclidean distance.
     * @throws TSException if error occures.
     */
    public static double seriesDistance(double[] series1, double[] series2) {
        assert series1.length == series2.length : "Exception in Euclidean distance: array lengths are not equal";
        Double res = 0D;
        for (int i = 0; i < series1.length; i++) {
            res = res + distance2(series1[i], series2[i]);
        }
        return Math.sqrt(res);
    }

    /**
     * Calculates euclidean distance between two multi-dimensional time-series
     * of equal length.
     *
     * @param series1 The first series.
     * @param series2 The second series.
     * @return The eclidean distance.
     * @throws TSException if error occures.
     */
    public static double seriesDistance(double[][] series1, double[][] series2) {
        assert series1.length == series2.length : "Exception in Euclidean distance: array lengths are not equal";
        Double res = 0D;
        for (int i = 0; i < series1.length; i++) {
            res = res + distance2(series1[i], series2[i]);
        }
        return Math.sqrt(res);

    }

    @Override
    public void clearCount() {
        cnt = 0;
    }

    @Override
    public long getCount() {
        return cnt;
    }
}
