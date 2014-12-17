PoreView
========

An open source Matlab data viewer based loosely on pClamp's ClampFit utility, 
along with a cached file loader for greatly simplified file access.



#Installation

Extract and copy to a folder, or fork the project via git. It is recommended that you add the folder to your Matlab
path to avoid any issues.

Tested in:
* Windows XP
* Windows 7
* Mac OSX




#Quick Start

PoreView is the GUI program that lets you view and manipulate large datafiles efficiently, as well as add dynamic
filters and superimpose plot objects.

## Launching

Recommended: `pv = pv_launch();` or `pv=pv_launch(filename);` or `pv=pv_launch(start_dir);`
Also works: `PoreView`, `pv=PoreView()`, `pv=PoreView(filename)` etc.

Think of the `pv_launch` file as a startup script - it sets up user-defined keyboard commands for the PoreView program, as well as 
starts it in a directory of your choosing, so that File->Open will go to the right place. By writing `pv = pv_launch`, you are saving a 
handle to the PoreView instance for later use, eg. `pv.autoscaleY()`.

Press 'Escape' to quit at any time.


## Adjusting view

* Scrolling over a plot zooms the x-axis about the cursor position. 
* Scrolling over a y-axis zooms that y axis about the cursor position.
* Middle-click and drag pans the plot.
* Right-clicking brings up the signal context menu.
* Left-click and drag on axes to zoom.
* Press 'a' to autoscale all axes.


## Cursors
* Double-click to bring cursors.
* Click and drag cursors to move them.
* Press 'c' to show/hide cursors.


## Useful functions

The following is defined by default in `pv_launch`:
* `f`: add filter easily
* `n`: plot noise, as in ClampFit
* `p`: plot visible signals in a new figure
* `s`: select raw data channels (for .fast5 raw data only)


## Plotting

For PoreView instance `pv`, you can do:

* `h = pv.getAxes(1)` returns handle to axes in topmost signal panel
* `plot(h, xdata, ydata)` plots on the upper panel
* `pv.clearAxes()` clears all user-plotted objects
* `pv.refresh()` forces redraw of all panels
* `pv.setView(trange)` moves visible window to `trange=[t0 t1]`


## SignalData

PoreView uses the SignalData class, which has a couple of features. The first time you load a file, it will build
a "reduced" version of the file, which contains the same data at a lower sampling rate, but using a type of min-max
sampling so that no spikes or peaks are lost. For large files, this can initially take a few minutes, but once saved they can be quickly loaded.
Another feature is caching - when you request points from SignalData, it only reloads from disk as necessary. It also doesn't crash when you request points past the end. This lets you do much more efficient sequential file processing.

The class can be accessed as `pv.data` for a PoreView instance `pv`.



## Data access

For PoreView instance `pv`, you can do:

* `trange = pv.getView()` to get visible data range
* `d = pv.data.getByTime(trange)` to get the data in that range
* `d = pv.data.getByTime(pv.getCursors())` to get the data in between the cursors
* `d = pv.data.get(range, cols)` returns data by indices in range, in specified columns only
* `d = pv.data.getViewData(trange)` returns reduced or full data based on zoom level.

For all data returned by SignalData, the first column is time and the rest are the other signals, real and filtered.

**Note:** `pv.data.get` and `pv.data.getByTime` can easily request enough data to make your computer unhappy! Use `getViewData` if you just want to look at the data.



## Data processing

For PoreView instance `pv`, you can do:

* `ind = pv.data.findNext(fun, istart)` finds the next time function 'fun' is true, from istart
* `ind = pv.data.findPrev(fun, istart)` reverse of above
* `pv.data.addVirtualSignal(fun, name)` adds a filter ('virtual signal') to the data

Examples:

* `ind = pv.data.findNext(@(d) d(:,2) > 0.5, 1000)` finds the first time the first signal (d(:,1) is time) is larger than 0.5, starting at index 1000
* `pv.data.addVirtualSignal(@(d) filt_lp(d,4,10000),'Low-pass 10kHz')` adds a 10 kHz low-pass filter to the data

Also refer to `find_events.m` to see an example of many of the SignalData and PoreView functions being used.


## More documentation

See the help for any of the classes PoreView and SignalData, or any of the functions, for more information. In Matlab, you can type `doc PoreView` to start the interactive help viewer at any time, or from the Help menu in PoreView.
