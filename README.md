# Perspective

## What
A script that finds ways to do perspective in .ass subtitles if you want rectangle to be rotated into a known tetragon.

## When
When you need to do a rotated sign and you can see or guess where rectangle should appear on the picture after rotation/perspective.

## Why
Because doing it by hand can be bothersome and because humans are frequently bad at this.

## How
0. Get python2, (python3 works if you replace "xrange" with "range")
1. Draw a tetragon that used to be a rectangle with Aegisub's clipping tool or get coordinates of four points otherwise. Points must be listed in clockwise order starting from top-left (relative to plane before rotation).
2. Get approximate \pos coordinates of your sign (or multiple coordinates for multiple signs) if you want to avoid using \org tag and try to guess width/height ratio of a unrotated rectangle (ratio of top and left sides of tetragon are used by default)
3. Call with coordinates of tetragon, optionally provide it with origin, width/height ratio and scaling factor (e.g. 0.66 if you have coordinates from 1920x1080 picture and want to typeset on 1280x720 picture). All arguments should be whatever-separated list of decimals, e.g. "{\clip(m 488 60 l 982, 360.7 990 1014 464 990)" or "\org(122.5,33);\org(422,133)"
4. Use some of resulting tags

## Which
For each specified target origin points script returns a perspective with \fax (when it's needed).
You're supposed to use it without \org and it will look fine as long as your sign is positioned near origin you passed to script.

Additionally script returns some proper perspectives with \org and without \fax.
Among possible origins it picks the one nearest to tetragon and ones that result in width/height ratio closer to target.
Former is usually easier to use because of closer origin.
This kind of perspective is easier to use for multiple signs or multiliners but it may be harder to position or use with motion tracking compared to perspective without \org. Don't try to use origin with coordinates over 100000000, renderers won't work with it properly.

### Hopefully easier-to-get instructions
![alt text](http://puu.sh/umBNo/d9c55343fa.png "old version of script, though")

### Extra thoughts: scaling rotation tags
You can probably project a rotation with simple formulas, rescale and unproject it with this script to achieve scaling of \frx and friends.
Alternatively, you can use something like this: [https://gist.github.com/Zeght/e997e663ca47d6aca2e7ee00dc92d62b](https://gist.github.com/Zeght/e997e663ca47d6aca2e7ee00dc92d62b)
