set -o xtrace
ORDER=(RandomOrientations Spins SpinsJoin SpinsJoinMain SpinRFPulse SpinRFPulseCoil)
parallel=5
for i in 1 2
do
    scenes=("${ORDER[@]:i:parallel}")
    echo $scenes
    for scene in ${scenes[@]}
    do
        manim -qh --renderer=opengl --write_to_movie spins.py $scene & 
    done
    wait
done
manim-slides convert --to=pptx ${ORDER[@]}
