reset
set cbrange [-10:10]
#set palette rgb 33,13,10
set palette model HSV
set palette rgb 3,2,2
set xrange [-0.5:99.5]
set yrange [-0.5:99.5]
unset xtics
unset ytics
do for [i=0:1000]{
    p "data.dat" i i matrix with image
    pause 0.05
}
