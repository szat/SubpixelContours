# SubpixelContours
A fast algorithm to create contours with a specified number of points or a specified mean distance (between points). Does not get confused when the contour has appendices of only one pixel thick. 

## Getting Started

This little project provides a resilient method for constructing contours with a subpixel precision. The user can pick the level of smoothness, and the total number of contour points or the average mean distance between the points. In particular, it does not get confused with contours that at some places are only one pixel thick. The caveat, is that it must be fed with 4-connected contours and not 8-connected contours.  

### Prerequisites

For this to run you need Matlab installed, along these toolboxes:
	-image_toolbox
	-statistics_toolbox

## Author

Adrian Szatmari 

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details