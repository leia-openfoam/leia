#!/usr/bin/awk -f
#
# Print each record EXCEPT
# IF the first record contains "raspberry", 
# THEN replace "red" with "pi"

BEGIN {
#	FS=",";
    OFS=","
    print "time", "symError"
}

/^Time =/{
    t=$3;
#    print t;
}
/symmetry-error/{
    se=$2;
    print t, se;
}
