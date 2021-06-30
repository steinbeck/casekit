package casekit.nmr.utils;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;

public class Statistics {
    /**
     * Detects outliers in given array list of input values and removes them. <br>
     * Here, outliers are those which are outside of a calculated lower and upper bound (whisker).
     * The interquartile range (IQR) of the input values is therefore multiplied with a given value
     * for whisker creation.
     *
     * @param input         list of values to process
     * @param multiplierIQR multiplier for IQR to use for lower and upper bound creation
     *
     * @return new array list without values outside the generated boundaries
     */
    public static List<Double> removeOutliers(final List<Double> input, final double multiplierIQR) {
        final ArrayList<Double> inputWithoutOutliers = new ArrayList<>(input);
        inputWithoutOutliers.removeAll(getOutliers(inputWithoutOutliers, multiplierIQR));

        return inputWithoutOutliers;
    }

    /**
     * @param input
     * @param multiplierIQR
     *
     * @return
     */
    public static List<Double> getOutliers(final List<Double> input, final double multiplierIQR) {
        final ArrayList<Double> outliers = new ArrayList<>();
        if (input.size()
                <= 1) {
            return outliers;
        }
        Collections.sort(input);
        final ArrayList<Double> data1 = new ArrayList<>(input.subList(0, input.size()
                / 2));
        final ArrayList<Double> data2;
        if (input.size()
                % 2
                == 0) {
            data2 = new ArrayList<>(input.subList(input.size()
                                                          / 2, input.size()));
        } else {
            data2 = new ArrayList<>(input.subList(input.size()
                                                          / 2
                                                          + 1, input.size()));
        }
        final double q1 = getMedian(new ArrayList<>(data1));
        final double q3 = getMedian(new ArrayList<>(data2));
        final double iqr = q3
                - q1;
        final double lowerBound = q1
                - multiplierIQR
                * iqr;
        final double upperBound = q3
                + multiplierIQR
                * iqr;
        for (int i = 0; i
                < input.size(); i++) {
            if ((input.get(i)
                    < lowerBound)
                    || (input.get(i)
                    > upperBound)) {
                outliers.add(input.get(i));
            }
        }
        //        System.out.println("input size: " + input.size());
        //        System.out.println("output size: " + outliers.size());
        return outliers;
    }

    /**
     * @param data
     *
     * @return
     */
    public static Double getMedian(final List<Double> data) {
        if ((data
                == null)
                || data.isEmpty()) {
            return null;
        }
        if (data.size()
                == 1) {
            return data.get(0);
        }
        Collections.sort(data);
        if (data.size()
                % 2
                == 1) {
            return data.get(data.size()
                                    / 2);
        } else {
            return (data.get(data.size()
                                     / 2
                                     - 1)
                    + data.get(data.size()
                                       / 2))
                    / 2.0;
        }
    }

    /**
     * @param data
     *
     * @return
     */
    public static Double getMean(final Collection<Double> data) {
        if ((data
                == null)
                || data.isEmpty()) {
            return null;
        }
        double sum = 0;
        int nullCounter = 0;
        for (final Double d : data) {
            if (d
                    != null) {
                sum += d;
            } else {
                nullCounter++;
            }
        }
        return ((data.size()
                - nullCounter)
                != 0)
               ? (sum
                / (data.size()
                - nullCounter))
               : null;
    }

    /**
     * @param data
     *
     * @return
     */
    public static Double getStandardDeviation(final List<Double> data) {
        if ((data
                == null)
                || data.isEmpty()) {
            return null;
        }
        final Double variance = getVariance(data);

        return (variance
                != null)
               ? Math.sqrt(variance)
               : null;
    }

    /**
     * @param data
     *
     * @return
     */
    public static Double getVariance(final Collection<Double> data) {
        if ((data
                == null)
                || data.isEmpty()) {
            return null;
        }
        final int nullCounter = Collections.frequency(data, null);
        double quadrSum = 0.0;
        final Double mean = getMean(data);
        if (mean
                == null) {
            return null;
        }
        for (final Double d : data) {
            if (d
                    != null) {
                quadrSum += Math.pow(d
                                             - mean, 2);
            }
        }

        return ((data.size()
                - nullCounter)
                != 0)
               ? (quadrSum
                / (data.size()
                - nullCounter))
               : null;
    }

    /**
     * @param data
     *
     * @return
     */
    public static Double getMean(final Double[] data) {
        if ((data
                == null)
                || (data.length
                == 0)) {
            return null;
        }
        double sum = 0;
        int nullCounter = 0;
        for (final Double d : data) {
            if (d
                    != null) {
                sum += d;
            } else {
                nullCounter++;
            }
        }
        return ((data.length
                - nullCounter)
                != 0)
               ? (sum
                / (data.length
                - nullCounter))
               : null;
    }

    public static Double roundDouble(final Double value, final int decimalPlaces) {
        if (value
                == null) {
            return null;
        }
        final int decimalFactor = (int) (Math.pow(10, decimalPlaces));

        return (Math.round(value
                                   * decimalFactor)
                / (double) decimalFactor);
    }

    /**
     * @param data
     *
     * @return
     */
    public static Double getRMSD(final Double[] data) {
        if (data
                == null
                || data.length
                == 0) {
            return null;
        }
        if (data.length
                == 1) {
            return data[0];
        }
        int nullCounter = 0;
        double qSum = 0;
        for (final Double d : data) {
            if (d
                    != null) {
                qSum += d
                        * d;
            } else {
                nullCounter++;
            }
        }

        return ((data.length
                - nullCounter)
                != 0)
               ? Math.sqrt(qSum
                                   / (data.length
                - nullCounter))
               : null;
    }

    /**
     * Returns the average of all deviations within a given input array.
     *
     * @param deviations array of deviations
     *
     * @return
     */
    public static Double calculateAverageDeviation(final Double[] deviations) {
        // every signal has to have a match
        for (final Double deviation : deviations) {
            if (deviation
                    == null) {
                return null;
            }
        }

        return getMean(deviations);
    }

    /**
     * Returns the average of all deviations within a given input array.
     *
     * @param data array of deviations
     *
     * @return
     */
    public static Double calculateRMSD(final Double[] data) {
        // every signal has to have a match
        for (final Double value : data) {
            if (value
                    == null) {
                return null;
            }
        }

        return getRMSD(data);
    }
}
