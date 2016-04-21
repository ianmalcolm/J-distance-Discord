/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package SAXFactory;

import weka.core.Attribute;
import weka.core.Instances;

/**
 *
 * @author Pavel Senin
 */
public class SAXFactory {

    /**
     * Convert real-valued series into symbolic representation.
     *
     * @param vals Real valued timeseries.
     * @param windowSize The PAA window size.
     * @param cuts The cut values array used for SAX transform.
     * @return The symbolic representation of the given real time-series.
     * @throws TSException If error occurs.
     */
    public static char[] getSaxVals(double[] vals, int windowSize, double[] cuts) {
        char[] saxVals;
        if (windowSize == cuts.length + 1) {
            saxVals = TSUtils.ts2String(TSUtils.zNormalize(vals), cuts);
        } else {
            saxVals = TSUtils.ts2String(TSUtils.zNormalize(TSUtils.paa(vals, cuts.length + 1)), cuts);
        }
        return saxVals;
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

    /**
     * This function implements SAX MINDIST function which uses alphabet based
     * distance matrix.
     *
     * @param a The SAX string.
     * @param b The SAX string.
     * @param distanceMatrix The distance matrix to use.
     * @return distance between strings.
     * @throws TSException If error occurs.
     */
    public static double saxMinDist(char[] a, char[] b, double[][] distanceMatrix) throws Exception {
        if (a.length == b.length) {
            double dist = 0.0D;
            for (int i = 0; i < a.length; i++) {
                if (Character.isLetter(a[i]) && Character.isLetter(b[i])) {
                    int numA = Character.getNumericValue(a[i]) - 10;
                    int numB = Character.getNumericValue(b[i]) - 10;
                    if (numA > 19 || numA < 0 || numB > 19 || numB < 0) {
                        throw new Exception("The character index greater than 19 or less than 0!");
                    }
                    double localDist = distanceMatrix[numA][numB];
                    dist += localDist;
                } else {
                    throw new Exception("Non-literal character found!");
                }
            }
            return dist;
        } else {
            throw new Exception("Data arrays lengths are not equal!");
        }
    }

    /**
     * Generic method to convert the milliseconds into the elapsed time string.
     *
     * @param start Start timestamp.
     * @param finish End timestamp.
     * @return String representation of the elapsed time.
     */
    public static String timeToString(long start, long finish) {
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
