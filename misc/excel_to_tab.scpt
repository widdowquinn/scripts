# excel_to_tab.scpt
#
# This script takes as input an Excel workbook containing one or more
# worksheets. It creates a new directory with the same name as the workbook,
# with the appended string _extracted. This directory contains a set of
# tab-separated plaintext files, one per worksheet. Each file has the same
# name as the corresponding worksheet, with the extension .tab.
#
# Installation:
# 
# This script is written in AppleScript, and is only expected to work on OSX.
#
# Place the excel_to_tab.scpt into your ~/Library/Scripts directory (create
# this directory if it does not exist).
#
# Open AppleScript Editor (in /Applications/Utilities) and open the General
# Preferences. Check the Show Script menu in menu bar setting, and close
# AppleScript Editor. You should now see the script symbol in the top menu
# bar.
#
# Usage:
#
# Click on the script symbol in the menu bar, and select the excel_to_tab
# option (it will be in the lower section). This will open a file selection
# dialog box. Select the appropriate Excel file, and click Choose. The script
# will generate the output directory in the same location as the Excel file.
#
# (c) Leighton Pritchard
# Authors: Leighton Pritchard
#
# Contact:
# leighton.pritchard@hutton.ac.uk
#
# Leighton Pritchard,
# Information and Computing Sciences,
# James Hutton Institute,
# Errol Road,
# Invergowrie,
# Dundee,
# DD6 9LH,
# Scotland,
# UK
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Select Excel workbook via dialogue
set theWorkbookFile to choose file with prompt "Please select an Excel Workbook"

# Open Excel and get useful information
tell application "Microsoft Excel"
  open theWorkbookFile
  set workbookName to name of active workbook
  if workbookName ends with ".xls" then set workbookName to text 1 thru -5 of workbookName
  if workbookName ends with ".xlsx" then set workbookName to text 1 thru -6 of workbookName
end tell

# Create new folder for output
set outputDirectory to (theWorkbookFile as text) & "_extracted"
if outputDirectory ends with ":" then set outputDirectory to text 1 thru -2 of outputDirectory
do shell script "mkdir -p " & quoted form of POSIX path of outputDirectory

# Loop over worksheets and write out in tab-separated format
tell application "Microsoft Excel"
  set theSheets to worksheets of active workbook
  repeat with aSheet in theSheets
    set thisPath to outputDirectory & ":" & workbookName & "_" & name of aSheet & ".tab"
    save aSheet in thisPath as text Mac file format
  end repeat
  close active workbook without saving
end tell
