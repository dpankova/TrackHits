# TrackHits
D.Pankova 03/2019

This is a part of StartingTrackVeto algorithm by K.Jero.
This part finds compatible hits for a given track.
Compatible Hits are the hits that could have come from that Track.
They are found using probabily tables:
1) You take a source - a segment of the track and a hit DOM
2) You Look into the tables to read out mean PE vs time 
   (given time interval and number of bins)
3) Find Max of mean PE vs time graph
4) Find two instances of time where mean PE = pecent*Max (default percent = 0.01). 
   This is time window, for compatible hits
5) Check if Pulses in the Dom occur during that time window 
   (then they are compatible)
6) Repeat for every segment and hit DOM

Unlike STV TrackHits module iterates over hit DOMs only, 
not all the doms withtin R (min_cad_dist). So it's faster.

Usage: look at /TrackHits/resources/test/RunTrackHits.py
