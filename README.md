#Perspective
## What
A script that finds ways to do perspective in .ass subtitles if you want rectangle to be rotated into a known tetragon.
## When
When you need do a rotated sign and you can see or guess where rectangle should appear on the picture after rotation/perspective.
## Why
Because doing it by hand can be bothersome and because humans are frequently bad at this.
## How
0. Get python2, (python3 works if you replace "xrange" with "range")
1. Draw a tetragon that used to be a rectangle with Aegisub's clipping tool or get four coordinates otherwise. Points must be listed in clockwise order starting from top-left (relative to plane before rotation).
2. Get approximate \pos coordinates your sign if you want to avoid using \org tag and try to guess width/height ratio of a unrotated rectangle (ratio of top and left sides of tetragon are used by default)
3. Call with a coordinates of tetragon, optionally provide it with origin and width/height ratio. All arguments should be whatever-separated list of floats, e.g. "{\clip(m 488 60 l 982, 360.7 990 1014 464 990)" or "\org(122,33)"
4. Use some of resulting tags

## Which
For specified target origin script returns a perspective with \fax (when it's needed).
You're supposed to use it without \org and it will look fine as long as your sign is positioned near origin you passed to script.

Additionally script returns some proper perspectives with \org and without \fax.
Among possible origins it picks the one closest to tetragon and ones that result in width/height ratio closer to target.
Former is usually easier to use because of closer origin.
This kind of perspective is easier to use for multiple signs or multiliners but it may be harder to position or use with motion tracking compared to perspective without \org.

### Hopefully easier-to-get instructions
![alt text](http://puu.sh/umBNo/d9c55343fa.png "old version of script, though")