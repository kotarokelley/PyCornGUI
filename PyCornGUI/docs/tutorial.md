#Tutorial for PyCornGUI
* Load PyCornGUI by double clicking on the app icon or from a terminal. Depending on where PyCornGUI is launched from, all .res files will be loaded. In this case, there are no .res files.  

![screen1](screenshots/screen_1.tiff)

* Load .res files from a directory by launching the Open pop-up by clicking on the Open Folder button.

* Before we get to displaying and saving traces, lets set the output directory by launching the Save Options pop-up by clicking on Save Options button from the options menu located at the tool bar. I usually set this to the same directory where the .res files are. This pop up also allows you to change the file fomat of the output. 

![screen2](screenshots/screen_2.tiff)

* Lets display a 280nm trace from the directory examples. To do this, high light the 280nm button by clicking on it at the bottom of the app, followed by double clicking on a file name from the File browser at the right side. This will display the entire 280nm trace for this run. 
![screen3](screenshots/screen_3.tiff)

* To display only a part of the trace, designate a range at the bottom of the app. The min and max has to be separated by a dash(-) or an underscore(_). Each time a parameter is changed, you have to redisplay by double clicking on the file name at the right. 

![screen4](screenshots/screen_4.tiff)

* To change display parameters, launch the Display Options pop-up from the tool bar at the top. Line color, thickness, legend, linear scaling(for multiple files), etc. 

![screen5](screenshots/screen_5.tiff)

* Once everything looks good, output trace to file by clicking Output selected to file button at the botom right. 

* If you have multiple files to process in the same way, then select multiple file names by clicking on them while holding down the command button on the keyboard, followed by clicking Output selected to file. As each file is processed, the trace will be displayed to the screen sequentially. 

* To overlay multiple traces, first select the file names by clicking on them while holding down the command button, follwed by clicking on Overlay selected to file button at the bottom right. Currently there is no option to view without saving so the overlay trace will be automatically saved. 

![screen6](screenshots/screen_6.tiff)

* As you can see, the two traces are displaced by a significant amount in the y axis. To mitigate this, select linear scaling check box from the Display Options pop-up and then re-displaying. 

![screen7](screenshots/screen_7.tiff)

* Thats about it. If a file contains a 260nm, 220nm, or conductance trace, they can be displayed by high lighting their corresponding buttons at the bottom. For some models of AKTA purifiers, the .res file may not contain any wavelength information. In those cases, display UV traces by high lighting the OtherUV button. This is also the case for other wavelengths such as for GFP. 
* Any requests for more features or if you notice something doesn't work as expected, would be greatly appreciated!


