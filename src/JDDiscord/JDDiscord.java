/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package JDDiscord;

import com.infomatiq.jsi.Point;
import com.infomatiq.jsi.Rectangle;
import com.infomatiq.jsi.SpatialIndex;
import com.infomatiq.jsi.detect.LOF;
import com.infomatiq.jsi.rtree.RTree;
import static edu.hawaii.jmotif.sax.SAXFactory.DEFAULT_COLLECTION_SIZE;
import static edu.hawaii.jmotif.sax.SAXFactory.getSaxVals;
import static edu.hawaii.jmotif.sax.SAXFactory.saxMinDist;
import edu.hawaii.jmotif.sax.alphabet.NormalAlphabet;
import edu.hawaii.jmotif.sax.datastructures.DiscordRecord;
import edu.hawaii.jmotif.sax.datastructures.DiscordRecords;
import edu.hawaii.jmotif.sax.trie.SAXTrie;
import edu.hawaii.jmotif.sax.trie.SAXTrieHitEntry;
import edu.hawaii.jmotif.sax.trie.TrieException;
import edu.hawaii.jmotif.timeseries.TSException;
import edu.hawaii.jmotif.timeseries.TSUtils;
import edu.hawaii.jmotif.util.BriefFormatter;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collections;
import java.util.Date;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.ListIterator;
import java.util.PriorityQueue;
import java.util.logging.ConsoleHandler;
import java.util.logging.Formatter;
import java.util.logging.Handler;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.hackystat.utilities.logger.HackystatLogger;
import weka.core.Attribute;
import weka.core.Instances;

/**
 *
 * @author ian
 */
public class JDDiscord {

    private static Logger consoleLogger;
    private static final String LOGGING_LEVEL = "SEVERE";

    static {
        consoleLogger = HackystatLogger.getLogger("jmotif.debug.console", "jmotif");
        consoleLogger.setUseParentHandlers(false);
        for (Handler handler : consoleLogger.getHandlers()) {
            consoleLogger.removeHandler(handler);
        }
        ConsoleHandler handler = new ConsoleHandler();
        Formatter formatter = new BriefFormatter();
        handler.setFormatter(formatter);
        consoleLogger.addHandler(handler);
        HackystatLogger.setLoggingLevel(consoleLogger, LOGGING_LEVEL);
    }

    public static DiscordRecords findBFDiscords(double[] series, int windowSize, int reportNum) throws TSException, TrieException {
        DiscordRecords discords = new DiscordRecords(reportNum);
        BitSet visitedDiscords = new BitSet(series.length - windowSize + 1);
        while (discords.getSize() < reportNum) {
            consoleLogger.log(Level.FINE, "currently known discords: {0} out of {1}", new Object[]{discords.getSize(), reportNum});
            Date start = new Date();
            DiscordRecord bestDiscord = getBFDiscord(series, windowSize, visitedDiscords);
            Date end = new Date();
            // if the discord is null we getting out of the search
            if (bestDiscord.getDistance() == 0.0D || bestDiscord.getPosition() == -1) {
                consoleLogger.log(Level.FINE, "breaking the outer search loop, discords found: {0} last seen discord: {1}", new Object[]{discords.getSize(), bestDiscord.toString()});
                break;
            }
//            consoleLogger.fine("new discord: " + bestDiscord.getPayload() + ", position "
//                    + bestDiscord.getPosition() + ", distance " + bestDiscord.getDistance()
//                    + ", elapsed time: " + timeToString(start.getTime(), end.getTime()));
            System.out.println("new discord: " + bestDiscord.getPayload() + ", position "
                    + bestDiscord.getPosition() + ", distance " + bestDiscord.getDistance()
                    + ", elapsed time: " + timeToString(start.getTime(), end.getTime()));
            // collect the result
            //
            {
                int pbegin = bestDiscord.getPosition() - windowSize + 1;
                int pend = bestDiscord.getPosition() + windowSize - 1;
                if (pbegin < 0) {
                    pbegin = 0;
                }
                if (pend >= visitedDiscords.size()) {
                    pend = visitedDiscords.size();
                }
                visitedDiscords.set(pbegin, pend);
            }
            discords.add(bestDiscord);
        }
        return discords;
    }

    /**
     *
     * "We are given n, the length of the discords in advance, and we must
     * choose two parameters, the cardinality of the SAX alphabet size a, and
     * the SAX word size w. We defer a discussion of how to set these parameters
     * until"
     *
     *
     * @param tsData timeseries.
     * @param windowLength window length.
     * @param visitedDiscords
     * @return top discords for the time-series given
     * @throws TSException if error occurs.
     */
    public static DiscordRecord getBFDiscord(double[] tsData, int windowLength, BitSet visitedDiscords)
            throws TSException {

        int bsfPosition = -1;
        double bsfDistance = Double.MIN_VALUE;

        // run the search loop
        //
        for (int i = 0; i < tsData.length - windowLength; i++) {
            if (visitedDiscords.get(i)) {
                continue;
            }

            double[] seriesA = getSubSeries(tsData, i, i + windowLength);
            Double nearestNeighborDist = Double.MAX_VALUE;

            // the inner loop
            //
            for (int j = 0; j < tsData.length - windowLength; j++) {

                if (Math.abs(i - j) < windowLength) {
                    continue;
                }

                double[] seriesB = getSubSeries(tsData, j, j + windowLength);
//                double dist = EuclideanDistance.distance(seriesA, seriesB);
                double dist = ZNDist2(seriesA, seriesB);

                if (dist < nearestNeighborDist) {
                    nearestNeighborDist = dist;
                }

            } // inner loop

            if (nearestNeighborDist > bsfDistance) {
                bsfDistance = nearestNeighborDist;
                bsfPosition = i;
            }
        }

        return new DiscordRecord(bsfPosition, bsfDistance);
    }

    public static DiscordRecords instanceBFDiscords(Instances tsData, String dataAttributeName,
            int windowSize, int discordsNumToReport) throws TrieException, TSException {

        Attribute dataAttribute = tsData.attribute(dataAttributeName);
        double[] series = toRealSeries(tsData, dataAttribute);

        Date start = new Date();
        int reportNum = DEFAULT_COLLECTION_SIZE;
        if (discordsNumToReport > 0 && discordsNumToReport < 50) {
            reportNum = discordsNumToReport;
        }
        DiscordRecords discords = findBFDiscords(series, windowSize, reportNum);
        Date end = new Date();
        consoleLogger.log(Level.FINE, "discords search finished in : {0}", timeToString(start.getTime(), end.getTime()));

        return discords;
    }

    public static DiscordRecords instanceJDDiscords(Instances tsData, String dataAttributeName,
            int windowSize, int alphabetSize, int discordsNumToReport, int j) throws TrieException, TSException {

        Attribute dataAttribute = tsData.attribute(dataAttributeName);
        double[] series = toRealSeries(tsData, dataAttribute);
        NormalAlphabet normalA = new NormalAlphabet();
        SAXTrieWithBound trie = new SAXTrieWithBound(series.length - windowSize, alphabetSize);

        StringBuilder sb = new StringBuilder();
        sb.append("data size: ").append(series.length);

        double max = TSUtils.max(series);
        sb.append("; max: ").append(max);

        double min = TSUtils.min(series);
        sb.append("; min: ").append(min);

        double mean = TSUtils.mean(series);
        sb.append("; mean: ").append(mean);

        int nans = TSUtils.countNaN(series);
        sb.append("; NaNs: ").append(nans);

        consoleLogger.fine(sb.toString());
        consoleLogger.log(Level.FINE, "window size: {0}, alphabet size: {1}, SAX Trie size: {2}", new Object[]{windowSize, alphabetSize, series.length - windowSize});

        // build the trie
        //
        Date start = new Date();
        int currPosition = 0;
        KNearestNeighbor[] knn = new KNearestNeighbor[series.length - windowSize];
        while ((currPosition + windowSize) < series.length) {
            // get the window SAX representation
            double[] subSeries = getSubSeries(series, currPosition, currPosition + windowSize);
            char[] saxVals = getSaxVals(subSeries, windowSize, normalA.getCuts(alphabetSize));
            // add result to the structure

            trie.put(String.valueOf(saxVals), currPosition, TSUtils.min(subSeries), TSUtils.max(subSeries));
            knn[currPosition] = new KNearestNeighbor(j, windowSize);
            // increment the position
            currPosition++;
        }

        Date end = new Date();
        System.out.println("Build Trie: , elapsed time: " + timeToString(start.getTime(), end.getTime()));

        start = new Date();
        int reportNum = DEFAULT_COLLECTION_SIZE;
        if (discordsNumToReport > 0 && discordsNumToReport < 50) {
            reportNum = discordsNumToReport;
        }
        DiscordRecords discords = findJDDiscords(toRealSeries(tsData, dataAttribute), windowSize, reportNum, j, trie, knn);
        end = new Date();
        consoleLogger.log(Level.FINE, "discords search finished in : {0}", timeToString(start.getTime(), end.getTime()));

        return discords;
    }

    public static DiscordRecords findJDDiscords(double[] series, int windowSize, int reportNum, int j, SAXTrieWithBound trie, KNearestNeighbor[] knn) throws TSException, TrieException {
        DiscordRecords discords = new DiscordRecords(reportNum);
        BitSet visitedDiscords = new BitSet(series.length - windowSize + 1);
        while (discords.getSize() < reportNum) {
            consoleLogger.log(Level.FINE, "currently known discords: {0} out of {1}", new Object[]{discords.getSize(), reportNum});

            Date start = new Date();
            DiscordRecord bestDiscord = findBestJDDiscord(series, visitedDiscords, windowSize, j, trie, knn);
            Date end = new Date();

            // if the discord is null we getting out of the search
            if (bestDiscord.getDistance() == 0.0D || bestDiscord.getPosition() == -1) {
                consoleLogger.log(Level.FINE, "breaking the outer search loop, discords found: {0} last seen discord: {1}", new Object[]{discords.getSize(), bestDiscord.toString()});
                break;
            }

//            consoleLogger.fine("new discord: " + bestDiscord.getPayload() + ", position "
//                    + bestDiscord.getPosition() + ", distance " + bestDiscord.getDistance()
//                    + ", elapsed time: " + timeToString(start.getTime(), end.getTime()));
            System.out.println("new discord: " + bestDiscord.getPayload() + ", position "
                    + bestDiscord.getPosition() + ", distance " + bestDiscord.getDistance()
                    + ", elapsed time: " + timeToString(start.getTime(), end.getTime()));
            // collect the result
            //
//            {
//                int pbegin = bestDiscord.getPosition() - windowSize + 1;
//                int pend = bestDiscord.getPosition() + windowSize - 1;
//                if (pbegin < 0) {
//                    pbegin = 0;
//                }
//                if (pend >= visitedDiscords.size()) {
//                    pend = visitedDiscords.size();
//                }
//                visitedDiscords.set(pbegin, pend);
//            }
            visitedDiscords.set(bestDiscord.getPosition());
            discords.add(bestDiscord);

        }
        return discords;
    }

    public static DiscordRecords instanceLOFDiscords(Instances tsData, String dataAttributeName,
            int windowSize, int discordsNumToReport, int j) throws TrieException, TSException {
        Attribute dataAttribute = tsData.attribute(dataAttributeName);
        double[] series = toRealSeries(tsData, dataAttribute);
        ArrayList<Rectangle> rects = new ArrayList<>();
        SpatialIndex si = new RTree();

        Date start = new Date();
        int currPosition = 0;
        while ((currPosition + windowSize) < series.length) {
            // get the subsequence as rectangle
            double[] subSeries = getSubSeries(series, currPosition, currPosition + windowSize);
            Point p = new Point();
            for (int i = 0; i < subSeries.length; i++) {
                p.add(subSeries[i]);
            }
            Rectangle r = new Rectangle(p, p);

            // add result to the index
            rects.add(r);
            si.add(r);
            currPosition++;
        }

        Date end = new Date();
        System.out.println("Build Trie: , elapsed time: " + timeToString(start.getTime(), end.getTime()));

        start = new Date();
        int reportNum = DEFAULT_COLLECTION_SIZE;
        if (discordsNumToReport > 0 && discordsNumToReport < 50) {
            reportNum = discordsNumToReport;
        }
        DiscordRecords discords = findLOFDiscords(rects, windowSize, reportNum, j, si);
        end = new Date();
        consoleLogger.log(Level.FINE, "discords search finished in : {0}", timeToString(start.getTime(), end.getTime()));

        return discords;
    }

    public static DiscordRecords findLOFDiscords(ArrayList<Rectangle> rects, int windowSize, int reportNum, int j, SpatialIndex si) throws TSException, TrieException {
        DiscordRecords discords = new DiscordRecords(reportNum);
        BitSet visitedDiscords = new BitSet(rects.size());
        while (discords.getSize() < reportNum) {
            consoleLogger.log(Level.FINE, "currently known discords: {0} out of {1}", new Object[]{discords.getSize(), reportNum});

            Date start = new Date();
            DiscordRecord bestDiscord = findBestLOFDiscord(rects, visitedDiscords, j, si);
            Date end = new Date();

            // if the discord is null we getting out of the search
            if (bestDiscord.getDistance() == 0.0D || bestDiscord.getPosition() == -1) {
                consoleLogger.log(Level.FINE, "breaking the outer search loop, discords found: {0} last seen discord: {1}", new Object[]{discords.getSize(), bestDiscord.toString()});
                break;
            }

//            consoleLogger.fine("new discord: " + bestDiscord.getPayload() + ", position "
//                    + bestDiscord.getPosition() + ", distance " + bestDiscord.getDistance()
//                    + ", elapsed time: " + timeToString(start.getTime(), end.getTime()));
            System.out.println("new discord: " + bestDiscord.getPayload() + ", position "
                    + bestDiscord.getPosition() + ", distance " + bestDiscord.getDistance()
                    + ", elapsed time: " + timeToString(start.getTime(), end.getTime()));
            // collect the result
            //
            {
                int pbegin = bestDiscord.getPosition() - windowSize + 1;
                int pend = bestDiscord.getPosition() + windowSize - 1;
                if (pbegin < 0) {
                    pbegin = 0;
                }
                if (pend >= visitedDiscords.size()) {
                    pend = visitedDiscords.size();
                }
                visitedDiscords.set(pbegin, pend);
                for (int i = pbegin; i < pend; i++) {
                    si.delete(rects.get(i));
                }
            }
            discords.add(bestDiscord);

        }
        return discords;

    }

    public static DiscordRecord findBestLOFDiscord(ArrayList<Rectangle> rects, BitSet visitedDiscords, int j, SpatialIndex si)
            throws TSException, TrieException {
        double bsfdist = 0;
        int bsfpos = 0;
        for (int i = 0; i < rects.size(); i++) {
            if (visitedDiscords.get(i)) {
                continue;
            }
            LOF lof = LOF.lof(rects.get(i).copys(), j, rects, si);
            if (bsfdist < lof.getfactor()) {
                bsfdist = lof.getfactor();
                bsfpos = i;
            }
        }

        return new DiscordRecord(bsfpos, bsfdist);
    }

    /**
     * Converts Instances into double array.
     *
     * @param tsData The instances data.
     * @param dataAttribute The attribute to use in conversion.
     * @return real-valued array.
     */
    public static double[] toRealSeries(Instances tsData, Attribute dataAttribute) {
        double[] vals = new double[tsData.numInstances()];
        for (int i = 0; i < tsData.numInstances(); i++) {
            vals[i] = tsData.instance(i).value(dataAttribute.index());
        }
        return vals;
    }

    private static double ZNDist2(double[] point1, double[] point2) throws TSException {
        if (point1.length == point2.length) {
            double p1mean = TSUtils.mean(point1);
            double p2mean = TSUtils.mean(point2);
            double p1stdv = TSUtils.stDev(point1);
            double p2stdv = TSUtils.stDev(point2);
            if (p1stdv == 0) {
                p1stdv = Double.MIN_VALUE;
            }
            if (p2stdv == 0) {
                p2stdv = Double.MIN_VALUE;
            }
//            p1stdv = 1.0;
//            p2stdv = 1.0;
            Double sum = 0D;
            for (int i = 0; i < point1.length; i++) {
//                double znp1 = (point1[i] - p1mean) / p1stdv;
//                double znp2 = (point2[i] - p2mean) / p2stdv;
                double znp1 = (point1[i] - p1mean);
                double znp2 = (point2[i] - p2mean);

                sum = sum + (znp2 - znp1) * (znp2 - znp1);
            }
            return Math.sqrt(sum);
        } else {
            throw new TSException("Exception in Euclidean distance: array lengths are not equal");
        }
    }

    /**
     *
     * "We are given n, the length of the discords in advance, and we must
     * choose two parameters, the cardinality of the SAX alphabet size a, and
     * the SAX word size w. We defer a discussion of how to set these parameters
     * until"
     *
     *
     * @param tsData time series.
     * @param reportedDiscords
     * @param windowLength window length.
     * @param j
     * @param trie
     * @return top discords for the time-series given
     * @throws TSException if error occurs.
     * @throws edu.hawaii.jmotif.sax.trie.TrieException
     */
    public static DiscordRecord findBestJDDiscord(double[] tsData, BitSet reportedDiscords, int windowLength, int j, SAXTrieWithBound trie, KNearestNeighbor[] knn)
            throws TSException, TrieException {

        Neighbor bsf = new Neighbor(-1, 0);
        int wordcount = 0;

        // run the search loop
        //
        ArrayList<SAXTrieHitEntry> frequencies = trie.getFrequencies();
        Collections.sort(frequencies);
        ArrayList<String> allWords = trie.getAllWords();
        BitSet visitedp = new BitSet(reportedDiscords.size());
        visitedp.clear();
        // initialze  visited p
        for (int i = 0; i < reportedDiscords.size(); i++) {
            if (reportedDiscords.get(i)) {
                int pl = i - windowLength + 1 >= 0 ? i - windowLength + 1 : 0;
                int pr = i + windowLength - 1 <= reportedDiscords.size() ? i + windowLength - 1 : reportedDiscords.size();
                for (int k = pl; k < pr; k++) {
                    visitedp.set(k);
                }
            }
        }

        // outter loop
        for (String word : allWords) {
            List<Integer> currentOccurences = trie.getOccurences(word.toCharArray());
//            System.out.printf("%d\n", currentOccurences.size());

            for (Integer p : currentOccurences) {

                // determine if p should be visited
                knn[p].validate(reportedDiscords);
                if (knn[p].getKDistance() * 2 < bsf.getDistance()) {
                    visitedp.set(p);
                    LinkedList<Neighbor> neighbors = knn[p].getNeighbors();
                    for (Neighbor n : neighbors) {
                        visitedp.set(n.getPosition());
                    }
                }
                if (knn[p].getKDistance() < bsf.getDistance() || knn[p].isComplete() || visitedp.get(p)) {
                    continue;
                }
//                if (knn[p].getKDistance() < bsf.getDistance() || knn[p].isComplete()) {
//                    continue;
//                }

                // begin visit p
                visitedp.set(p);
                double[] seriesA = getSubSeries(tsData, p, p + windowLength - 1);
                BitSet visitedq = new BitSet(reportedDiscords.size());
                visitedq.clear();
                boolean completeSearch1 = true;
                boolean completeSearch2 = true;

                // inner loop
                // pick up the same pattern to compare
                for (Integer q : currentOccurences) {
                    if (Math.abs(q - p) < windowLength) {
                        continue;
                    }
                    // get the piece of the timeseries
                    visitedq.set(q);
                    double[] seriesB = getSubSeries(tsData, q, q + windowLength - 1);
//                    double dist = EuclideanDistance.distance(seriesA, seriesB);
                    double dist = ZNDist2(seriesA, seriesB);
                    wordcount++;
                    knn[p].add(q, dist);

                    if (knn[p].getKDistance() * 2 < bsf.getDistance()) {
                        LinkedList<Neighbor> neighbors = knn[p].getNeighbors();
                        for (Neighbor n : neighbors) {
                            visitedp.set(n.getPosition());
                        }
                        completeSearch1 = false;
                    }
                    if (knn[p].getKDistance() < bsf.getDistance()) {
                        completeSearch1 = false;
                        break;
                    }
                }

                if (completeSearch1 == true) {
                    // random pick the rest
                    ArrayList<Integer> sslist = trie.getAllSubsequences();
                    Collections.shuffle(sslist);
                    for (Integer q : sslist) {
                        if (visitedq.get(q)) {
                            continue;
                        }
                        if (Math.abs(q - p) < windowLength) {
                            continue;
                        }

                        // get the piece of the timeseries
                        visitedq.set(q);
                        double[] seriesB = getSubSeries(tsData, q, q + windowLength - 1);
//                    double dist = EuclideanDistance.distance(seriesA, seriesB);
                        double dist = ZNDist2(seriesA, seriesB);
                        wordcount++;
                        knn[p].add(q, dist);

                        if (knn[p].getKDistance() * 2 < bsf.getDistance()) {
                            LinkedList<Neighbor> neighbors = knn[p].getNeighbors();
                            for (Neighbor n : neighbors) {
                                visitedp.set(n.getPosition());
                            }
                            completeSearch2 = false;
                        }
                        if (knn[p].getKDistance() < bsf.getDistance()) {
                            completeSearch2 = false;
                            break;
                        }
                    }
                }

                // check if complete search J nearest neighbours of p
                if (completeSearch2) {
                    knn[p].setComplete();
                }

                // update best so far J-distance
                if (knn[p].getKDistance() > bsf.getDistance()) {
                    bsf = new Neighbor(p, knn[p].getKDistance());
                }
            }

        }

        System.out.println("Call distance function = " + wordcount);

        return new DiscordRecord(bsf.getPosition(), bsf.getDistance());
    }

    /**
     * Extracts sub-series from series.
     *
     * @param data The series.
     * @param start The start position.
     * @param end The end position
     * @return sub-series from start to end.
     */
    public static double[] getSubSeries(double[] data, int start, int end) {
        double[] vals = new double[end - start];
        for (int i = 0; i < end - start; i++) {
            vals[i] = data[start + i];
        }
        return vals;
    }

    /**
     * Generic method to convert the milliseconds into the elapsed time string.
     *
     * @param start Start timestamp.
     * @param finish End timestamp.
     * @return String representation of the elapsed time.
     */
    private static String timeToString(long start, long finish) {
        long diff = finish - start;

        long secondInMillis = 1000;
        long minuteInMillis = secondInMillis * 60;
        long hourInMillis = minuteInMillis * 60;
        long dayInMillis = hourInMillis * 24;
        long yearInMillis = dayInMillis * 365;

        @SuppressWarnings("unused")
        long elapsedYears = diff / yearInMillis;
        diff = diff % yearInMillis;

        @SuppressWarnings("unused")
        long elapsedDays = diff / dayInMillis;
        diff = diff % dayInMillis;

        @SuppressWarnings("unused")
        long elapsedHours = diff / hourInMillis;
        diff = diff % hourInMillis;

        long elapsedMinutes = diff / minuteInMillis;
        diff = diff % minuteInMillis;

        long elapsedSeconds = diff / secondInMillis;
        diff = diff % secondInMillis;

        long elapsedMilliseconds = diff % secondInMillis;

        return elapsedMinutes + "m " + elapsedSeconds + "s " + elapsedMilliseconds + "ms";

    }
}

class SAXTrieWithBound extends SAXTrie {

    private ArrayList<Double> seqmin = new ArrayList<>();
    private ArrayList<Double> seqmax = new ArrayList<>();
    private NormalAlphabet normalA = new NormalAlphabet();

    public SAXTrieWithBound(int seqLength, int alphabetSize) throws TrieException {
        super(seqLength, alphabetSize);
        for (int i = 0; i < seqLength; i++) {
            seqmin.add(Double.MAX_VALUE);
            seqmax.add(Double.MIN_VALUE);
        }
    }

    public PriorityQueue<Neighbor> genDistPQ(List<Integer> occu, double[] series) {
        PriorityQueue<Neighbor> distPQ = new PriorityQueue<>();
        double qmin = TSUtils.min(series);
        double qmax = TSUtils.max(series);
        for (int i = 0; i < occu.size(); i++) {
            int o = occu.get(i);
            Neighbor n = new Neighbor(o, qdist(seqmin.get(o), seqmax.get(o), qmin, qmax));
            distPQ.add(n);
        }
        return distPQ;
    }

    /**
     * Put an actual entrance in the magic trie.
     *
     * @param str The string.
     * @param idx The position.
     * @throws TrieException If error occurs.
     */
    public void put(String str, int idx, double min, double max) throws TrieException {
        super.put(str, idx);
        seqmin.set(idx, min);
        seqmax.set(idx, max);
    }

    public void sort(List<Integer> occu, double[] series) {
        double qmin = TSUtils.min(series);
        double qmax = TSUtils.max(series);
        ArrayList<Double> dist = new ArrayList();
        for (int i = 0; i < occu.size(); i++) {
            int o = occu.get(i);
            dist.add(qdist(seqmin.get(o), seqmax.get(o), qmin, qmax));
        }
        for (int i = 0; i < occu.size() - 1; i++) {
            for (int j = 0; j < occu.size() - i - 1; j++) {
                if (dist.get(j) > dist.get(j + 1)) {
                    int oswap = occu.get(j);
                    double dswap = dist.get(j);
                    occu.set(j, occu.get(j + 1));
                    dist.set(j, dist.get(j + 1));
                    occu.set(j + 1, oswap);
                    dist.set(j + 1, dswap);
                }
            }
        }
    }

    private double qdist(double ax, double ay, double bx, double by) {
        return (ax - bx) * (ax - bx) + (ay - by) * (ay - by);
    }

    public ArrayList<String> getAllWords() {
        ArrayList<SAXTrieHitEntry> freq = getFrequencies();
        ArrayList<String> wordList = new ArrayList();
        for (SAXTrieHitEntry a : freq) {
            String word = new String(a.getStr());
            if (wordList.contains(word)) {
                continue;
            }
            wordList.add(word);
        }
        return wordList;
    }

    public ArrayList<Integer> getAllSubsequences() {
        ArrayList<SAXTrieHitEntry> freq = getFrequencies();
        ArrayList<Integer> sslist = new ArrayList();
        for (SAXTrieHitEntry a : freq) {
            sslist.add(a.getPosition());
        }
        return sslist;
    }

    public void sortWords(ArrayList<String> wordList, String word) throws TSException {

        double[][] distMatrix = normalA.getDistanceMatrix(getAlphabetSize());
        ArrayList<Double> distList = new ArrayList();
        for (int i = 0; i < wordList.size(); i++) {
            distList.add(saxMinDist(wordList.get(i).toCharArray(), word.toCharArray(), distMatrix));
        }
        for (int i = 0; i < distList.size() - 1; i++) {
            for (int j = 0; j < distList.size() - i - 1; j++) {
                if (distList.get(j) > distList.get(j + 1)) {
                    String oswap = wordList.get(j);
                    double dswap = distList.get(j);
                    wordList.set(j, wordList.get(j + 1));
                    distList.set(j, distList.get(j + 1));
                    wordList.set(j + 1, oswap);
                    distList.set(j + 1, dswap);
                }
            }
        }
    }
}

class KNearestNeighbor implements Iterable<Neighbor> {

    private LinkedList<Neighbor> neighbors;
    private int maxCapacity;
    private boolean complete;
    private int windowLength;

    public KNearestNeighbor(int capacity, int length) {
        this.maxCapacity = capacity;
        this.neighbors = new LinkedList<>();
        this.complete = false;
        this.windowLength = length;
    }

    @Override
    public Iterator<Neighbor> iterator() {
        return this.neighbors.iterator();
    }

    public boolean isComplete() {
        return this.complete;
    }

//    public boolean isValid(BitSet vDs) {
//        for (Neighbor n : neighbors) {
//            if (vDs.get(n.getPosition())) {
//                return false;
//            }
//        }
//        return true;
//    }
    public void validate(BitSet vDs) {
        ListIterator<Neighbor> nli = neighbors.listIterator();
        while (nli.hasNext()) {
            if (vDs.get(nli.next().getPosition())) {
                nli.remove();
            }
        }
//        for (Neighbor n : neighbors) {
//            if (vDs.get(n.getPosition())) {
//                neighbors.remove(n);
//            }
//        }
    }

    public void setComplete() {
        this.complete = true;
    }

    public void clearComplete() {
        this.complete = false;
    }

    public boolean contains(int id) {
        for (int i = 0; i < neighbors.size(); i++) {
            if (neighbors.get(i).getPosition() == id) {
                return true;
            }
        }
        return false;
    }

    public int qualify(int id, double dist) {
        int flag = 0;
        int count = 0;

        for (int i = 0; i < neighbors.size(); i++) {
            if (Math.abs(neighbors.get(i).getPosition() - id) < this.windowLength) {
                count++;
                if (neighbors.get(i).getDistance() > dist) {
                    flag = i;
                }
            }
        }
        if (count > 1) {
            flag = -1;
        }
        return flag;
    }

    public boolean overlap(int id) {
        for (Neighbor n : neighbors) {
            if (Math.abs(n.getPosition() - id) < this.windowLength) {
                return true;
            }
        }
        return false;
    }

    public boolean add(int id, double dist) {
        if (!overlap(id)) {
            if (this.neighbors.size() < this.maxCapacity) {
                Neighbor n = new Neighbor(id, dist);
                this.neighbors.add(n);

            } else {
                if (neighbors.get(0).getDistance() > dist) {
                    Neighbor n = new Neighbor(id, dist);
                    neighbors.remove(0);
                    neighbors.add(n);

                }
            }
            Collections.sort(neighbors, Collections.reverseOrder());
            return true;
        } else {
            int pos = qualify(id, dist);
            if (pos > 0) {
                neighbors.remove(pos);
                neighbors.add(new Neighbor(id, dist));
            }
            Collections.sort(neighbors, Collections.reverseOrder());
            return false;
        }
    }

//    public boolean add(int id, double dist) {
//        if (!contains(id)) {
//            if (this.neighbors.size() < this.maxCapacity) {
//                Neighbor n = new Neighbor(id, dist);
//                this.neighbors.add(n);
//                return true;
//            } else {
//                if (neighbors.get(0).getDistance() > dist) {
//                    Neighbor n = new Neighbor(id, dist);
//                    neighbors.remove(0);
//                    neighbors.add(n);
//                    return true;
//                }
//            }
//            Collections.sort(neighbors, Collections.reverseOrder());
//        }
//        return false;
//    }
    public double getKDistance() {
        if (isFull()) {
            return neighbors.get(0).getDistance();
        } else {
            return Double.MAX_VALUE;
        }
    }

    public int getKPosition() {
        if (isFull()) {
            return neighbors.get(0).getPosition();
        } else {
            return -1;
        }
    }

    public Neighbor getKNeighbor() {

        if (isFull()) {
            return neighbors.get(0);
        } else {
            return null;
        }
    }

    public LinkedList<Neighbor> getNeighbors() {
        return this.neighbors;
    }

    public boolean isFull() {
        return neighbors.size() >= this.maxCapacity;
    }

}

class Neighbor implements Comparable<Neighbor> {

    private final int position;
    private final double distance;

    public Neighbor() {
        this.position = -1;
        this.distance = 0;
    }

    public Neighbor(int position, double distance) {
        this.position = position;
        this.distance = distance;
    }

    public double getDistance() {
        return this.distance;
    }

//    public void setDistance(double dist){
//        this.distance = dist;
//    }
    public int getPosition() {
        return this.position;
    }

    /**
     * The simple comparator based on the distance.
     *
     * @param other The discord record this is compared to.
     * @return True if equals.
     */
    @Override
    public int compareTo(Neighbor other) {
        if (null == other) {
            throw new NullPointerException("Unable compare to null!");
        }
        if (this.distance > other.getDistance()) {
            return 1;
        } else if (this.distance < other.getDistance()) {
            return -1;
        }
        return 0;
    }
}
