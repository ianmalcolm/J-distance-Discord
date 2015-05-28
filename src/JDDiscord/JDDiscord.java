/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package JDDiscord;

import SAXFactory.DiscordRecord;
import SAXFactory.DiscordRecords;
import SAXFactory.NormalAlphabet;
import SAXFactory.SAXFactory;
import SAXFactory.TSUtils;
import edu.hawaii.jmotif.sax.trie.SAXTrie;
import edu.hawaii.jmotif.sax.trie.SAXTrieHitEntry;
import edu.hawaii.jmotif.sax.trie.TrieException;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collections;
import java.util.Date;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.PriorityQueue;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author ian
 */
public class JDDiscord {

    private static final Logger jdlogger = Logger.getLogger(JDDiscord.class.getName());
    private final SAXTrieWithBound trie;
    private double[] series;
    private int windowSize;
    private int alphabetSize;
    private int dimension;
    private DataHandler dh;
    private Distance df;
    private boolean ENHANCED = true;
//    private ArrayList<String> outerWords;
    private ArrayList<String> innerWords;
    private ArrayList<SAXTrieHitEntry> frequencies;

    public long totalcnt = 0;

    public JDDiscord(double[] ts, int ws, int as, int dim, DataHandler _dh, Distance _df) throws TrieException {
        series = ts;
        windowSize = ws;
        alphabetSize = as;
        dimension = dim;
        dh = _dh;
        df = _df;

        NormalAlphabet normalA = new NormalAlphabet();
        trie = new SAXTrieWithBound(series.length - windowSize + 1, alphabetSize);

        StringBuilder sb = new StringBuilder();
        sb.append("data size: ").append(series.length);

        double max = TSUtils.max(series);
        sb.append("; max: ").append(max);

        double min = TSUtils.min(series);
        sb.append("; min: ").append(min);

        double mean = TSUtils.mean(series);
        sb.append("; mean: ").append(mean);

        double std = TSUtils.stDev(series);
        sb.append("; std: ").append(std);

        jdlogger.fine(sb.toString());
        jdlogger.fine("window size: " + windowSize + ", alphabet size: " + alphabetSize
                + ", SAX Trie size: " + (series.length - windowSize));

        // build the trie
        for (int i = 0; i < dh.size(); i++) {
            // get the window SAX representation
            char[] saxVals = SAXFactory.getSaxVals(dh.get(i), windowSize, normalA.getCuts(alphabetSize));
            // add result to the structure
            trie.put(String.valueOf(saxVals), i);
            // increment the position            
        }

//        outerWords = trie.getAllWords();
        innerWords = trie.getAllWords();
        frequencies = trie.getFrequencies();
        Collections.sort(frequencies);

    }

    public static void setLoggerLevel(Level level) {
        jdlogger.setLevel(level);
    }

    public void setEnhanced(boolean e) {
        ENHANCED = e;
    }

    public DiscordRecords findJDDiscords(int reportNum, int j) throws Exception {
        DiscordRecords discords = new DiscordRecords(reportNum);
        BitSet exception = new BitSet(series.length - windowSize + 1);
        KNearestNeighbor[] knn = new KNearestNeighbor[series.length - windowSize + 1];

        for (int i = 0; i < dh.size(); i++) {
            knn[i] = new KNearestNeighbor(j, windowSize);
        }

        while (discords.getSize() < reportNum) {
            jdlogger.fine("currently known discords: " + discords.getSize() + " out of " + reportNum);

            Date start = new Date();
            DiscordRecord bestDiscord = findBestJDDiscord(exception, windowSize, j, knn);
            Date end = new Date();

            // if the discord is null we getting out of the search
            if (bestDiscord.getDistance() < 0.0D || bestDiscord.getPosition() < 0) {
                jdlogger.fine("breaking the outer search loop, discords found: " + discords.getSize() + " last seen discord: " + bestDiscord.toString());
            }

            jdlogger.info(
                    "new discord: " + bestDiscord.getPayload()
                    + ", position " + bestDiscord.getPosition()
                    + ", distance " + bestDiscord.getDistance()
                    + ", elapsed time: " + SAXFactory.timeToString(start.getTime(), end.getTime())
                    + ", call to distance function: " + df.getCount());
            totalcnt += df.getCount();
            df.clearCount();
            // collect the result
            {
                int pbegin = (int) (bestDiscord.getPosition() - windowSize + 1);
                int pend = (int) (bestDiscord.getPosition() + windowSize - 1);
                if (pbegin < 0) {
                    pbegin = 0;
                }
                if (pend >= exception.size()) {
                    pend = exception.size();
                }
                exception.set(pbegin, pend);
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
     * @param exception
     * @param windowLength window length.
     * @param j
     * @param knn
     * @param trie
     * @return top discords for the time-series given
     * @throws Exception if error occurs.
     * @throws edu.hawaii.jmotif.sax.trie.TrieException
     */
    public DiscordRecord findBestJDDiscord(BitSet exception, int windowLength, int j, KNearestNeighbor[] knn) throws TrieException, Exception {

        Neighbor bsf = new Neighbor(Integer.MIN_VALUE, Double.NEGATIVE_INFINITY);

        // run the search loop
        BitSet visitedp = (BitSet) exception.clone();

        // outer loop
        for (SAXTrieHitEntry currentOccurences : frequencies) {
            int p = currentOccurences.getPosition();
            String outerWord = new String(currentOccurences.getStr());
            jdlogger.finer("Position: " + p + "\tload: " + outerWord + "\tfrequency: " + currentOccurences.getFrequency());
            // determine if p should be visited
            {
                if (ENHANCED) {
                    if (knn[p].getKDistance() * 2 < bsf.getDistance() && knn[p].getKDistance() > 0) {
                        visitedp.set(p);
                        knn[p].mark(visitedp);
                    }
                }
                if (knn[p].getKDistance() < bsf.getDistance() && knn[p].getKDistance() > 0) {
                    continue;
                }
                if (visitedp.get(p)) {
                    continue;
                }
            }

            // begin visit p
            visitedp.set(p);
            List<Integer> currentOuterOccurences = trie.getOccurences(outerWord.toCharArray());
            double[] seriesA = dh.get(p);
            boolean completeSearch = true;

            // inner loop
            for (Integer q : currentOuterOccurences) {
                if (Math.abs(q - p) < windowLength) {
                    continue;
                }

                // get the piece of the timeseries
                if (knn[p].contains(q)) {
                    continue;
                }
                double[] seriesB = dh.get(q);
                double dist = df.distance(seriesA, seriesB);
                knn[p].add(q, dist);
                knn[q].add(p, dist);

                {

                    if (ENHANCED) {
                        if (knn[p].getKDistance() * 2 < bsf.getDistance() && knn[p].getKDistance() > 0) {
                            visitedp.set(p);
                            knn[p].mark(visitedp);
                        }
                    }
                    if (knn[p].getKDistance() < bsf.getDistance() && knn[p].getKDistance() > 0) {
                        completeSearch = false;
                        break;
                    }
                }
            }

            if (completeSearch) {

                trie.sortWords(innerWords, outerWord);
                for (String innerWord : innerWords) {
                    if (innerWord.compareTo(outerWord) == 0) {
                        continue;
                    }
                    List<Integer> currentInnerOccurences = trie.getOccurences(innerWord.toCharArray());
                    for (Integer q : currentInnerOccurences) {
                        if (Math.abs(q - p) < windowLength) {
                            continue;
                        }

                        if (knn[p].contains(q)) {
                            continue;
                        }
                        double[] seriesB = dh.get(q);
                        double dist = df.distance(seriesA, seriesB);
                        knn[p].add(q, dist);
                        knn[q].add(p, dist);

                        {
                            if (ENHANCED) {
                                if (knn[p].getKDistance() * 2 < bsf.getDistance() && knn[p].getKDistance() > 0) {
                                    visitedp.set(p);
                                    knn[p].mark(visitedp);
                                }
                            }
                            if (knn[p].getKDistance() < bsf.getDistance() && knn[p].getKDistance() > 0) {
                                completeSearch = false;
                                break;
                            }
                        }
                    }

                    if (completeSearch == false) {
                        break;
                    }
                }
            }

            // update best so far J-distance
            if (knn[p].getKDistance() > bsf.getDistance()) {
                bsf = new Neighbor(p, knn[p].getKDistance());
                jdlogger.fine("update best so far position: " + bsf.getPosition() + "\tdistance: " + bsf.getDistance());
            }
        }

        return new DiscordRecord(bsf.getPosition(), bsf.getDistance());
    }

}

class SAXTrieWithBound extends SAXTrie {

    private ArrayList<Double> seqmin = new ArrayList<>();
    private ArrayList<Double> seqmax = new ArrayList<>();
    private NormalAlphabet normalA = new NormalAlphabet();
    private int size;

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
    public void put(String str, int idx, double min, double max) throws Exception {
        super.put(str, idx);
        seqmin.set(idx, min);
        seqmax.set(idx, max);
        size++;
    }

    public int size() {
        return size;
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

    public void sortWords(ArrayList<String> wordList, String word) throws Exception {

        double[][] distMatrix = normalA.getDistanceMatrix(getAlphabetSize());
        ArrayList<Double> distList = new ArrayList();
        for (int i = 0; i < wordList.size(); i++) {
            distList.add(SAXFactory.saxMinDist(wordList.get(i).toCharArray(), word.toCharArray(), distMatrix));
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
    private int windowSize;
    private double kdist = Double.NEGATIVE_INFINITY;

    public KNearestNeighbor(int capacity, int length) {
        this.maxCapacity = capacity;
        this.neighbors = new LinkedList();
        this.windowSize = length;
    }

    @Override
    public Iterator<Neighbor> iterator() {
        return this.neighbors.iterator();
    }

    private boolean contains(LinkedList<Neighbor> list, int id) {
        for (int i = 0; i < list.size(); i++) {
            if (list.get(i).getPosition() == id) {
                return true;
            }
        }
        return false;
    }

    public boolean contains(int id) {
        return contains(neighbors, id);
    }

    public boolean overlap(LinkedList<Neighbor> list, int id) {
        for (Neighbor n : list) {
            if (Math.abs(n.getPosition() - id) < this.windowSize) {
                return true;
            }
        }
        return false;
    }

    public boolean add(int id, double dist) {
        double _kdist = getKDistance();
        if (_kdist > dist || _kdist == Double.NEGATIVE_INFINITY) {
            Neighbor n = new Neighbor(id, dist);
            int position = Collections.binarySearch(neighbors, n);
            if (position < 0) {
                neighbors.add(-1 * (position + 1), n);
                return true;
            }
        }
        return false;
    }

    public double getKDistance() {

        if (neighbors.size() > 0) {
            LinkedList<Neighbor> knn = getKNN();
            if (knn.size() >= this.maxCapacity) {
                kdist = knn.get(knn.size() - 1).getDistance();
            }
        }

        return kdist;
    }

    public LinkedList<Neighbor> getKNN() {
        if (neighbors.size() > 0) {
            LinkedList<Neighbor> result = new LinkedList();
            for (Neighbor n : neighbors) {
                if (overlap(result, n.getPosition())) {
                    continue;
                }
                if (result.size() < maxCapacity) {
                    result.add(n);
                } else if (n.getDistance() == result.getLast().getDistance()) {
                    result.add(n);
                } else {
                    break;
                }
            }
            return result;
        } else {
            return null;
        }
    }

    public void mark(BitSet visited) {
        for (Neighbor n : neighbors) {
            visited.set(n.getPosition());
        }
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
