/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package JDDiscord;

/**
 *
 * @author ian
 */
public abstract class Distance {
    public abstract double distance(double[] ts1, double[] ts2);

    public abstract void clearCount();
    
    public abstract long getCount();
}
