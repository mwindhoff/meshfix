#!/bin/sh
../meshfix sphere1.off sphere2.off -n 2 -j 20 -jc -o join_two_overlapping_spheres_result
../meshfix close_objects.off -n 4 -j 30 -jc -o join_close_objects_result

