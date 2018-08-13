# Copyright (c) 2018, Gabriel Johnson. All rights reserved.

/^[[:space:]]*[[:digit:]]\{1,\}$/ {
N
s_\n[^-]*-\{1,\}_/_
N
s/\n//
s/[[:space:]]*//g
p
}
