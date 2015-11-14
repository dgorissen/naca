NACA [![Travis][buildstatus_image_travis]][travisci]

Code to generate 4 and 5 digit NACA profiles. 

Problems/Improvements: Dirk Gorissen <dgorissen@gmail.com>

    # Get help
    python naca.py -h
    # Generate points for NACA profile 2412
    python naca.py -p 2412
    # Generate points for NACA profile 2412 with 300 points
    python naca.py -p 2412 -n 300
    # Generate points for NACA profile 2412 and display the result
    python naca.py -p 2412 -d
    # Generate points for NACA profile 2412 with smooth points spacing and display the result
    python naca.py -p 2412 -d -s
    # Generate points for several profiles
    python naca.py -p "2412 23112" -d -s