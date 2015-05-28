/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package SAXFactory;


/**
 * The discord data record. Kind of useful of keeping information concerning a discord occurrence
 * etc...
 * 
 * @author Pavel Senin
 * 
 */
public class DiscordRecord implements Comparable<DiscordRecord> {

  /** Indicates the position of this discord on the timeline. */
  private long position;

  /** The distance of this discord instance to the closest neighbor. */
  private double distance;

  /** The payload. */
  private String payload;

  /**
   * Constructor.
   */
  public DiscordRecord() {
    this.position = -1;
    this.distance = 0D;
    this.payload = "";
  }

  /**
   * Constructor.
   * 
   * @param index The index discord found at.
   * @param dist The distance from other sequences.
   */
  public DiscordRecord(long index, double dist) {
    this.position = index;
    this.distance = dist;
    this.payload = "";
  }

  /**
   * Constructor.
   * 
   * @param index The index discord found at.
   * @param dist The distance from other sequences.
   * @param payload The payload.
   */
  public DiscordRecord(int index, double dist, String payload) {
    this.position = index;
    this.distance = dist;
    this.payload = payload;
  }

  /**
   * Set the payload value.
   * 
   * @param payload The payload.
   */
  public void setPayload(String payload) {
    this.payload = payload;
  }

  /**
   * Get the payload.
   * 
   * @return The payload.
   */
  public String getPayload() {
    return payload;
  }

  /**
   * Set the position at the time series list.
   * 
   * @param position the position to set.
   */
  public void setIndex(int position) {
    this.position = position;
  }

  /**
   * Get the position.
   * 
   * @return the position at the time series list.
   */
  public long getPosition() {
    return this.position;
  }

  /**
   * Set the distance to the closest neighbor.
   * 
   * @param distance the distance to set.
   */
  public void setDistance(double distance) {
    this.distance = distance;
  }

  /**
   * Get the distance to the closest neighbor.
   * 
   * @return the distance to the closest neighbor.
   */
  public double getDistance() {
    return this.distance;
  }

  /**
   * The simple comparator based on the distance.
   * 
   * @param other The discord record this is compared to.
   * @return True if equals.
   */
  @Override
  public int compareTo(DiscordRecord other) {
    if (null == other) {
      throw new NullPointerException("Unable compare to null!");
    }
    if (this.distance > other.getDistance()) {
      return 1;
    }
    else if (this.distance < other.getDistance()) {
      return -1;
    }
    return 0;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public boolean equals(Object obj) {
    if (this == obj) {
      return true;
    }
    if (obj == null) {
      return false;
    }
    if (getClass() != obj.getClass()) {
      return false;
    }
    DiscordRecord other = (DiscordRecord) obj;
    if (Double.doubleToLongBits(distance) != Double.doubleToLongBits(other.distance)) {
      return false;
    }
    if (position != other.position) {
      return false;
    }
    return true;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public String toString() {
    return "'" + this.payload + "', distance " + this.getDistance() + " position: "
        + this.getPosition();
  }

}
