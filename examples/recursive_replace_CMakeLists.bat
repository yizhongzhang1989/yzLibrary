@echo off

rem	this batch recursively replace all CMakeLists.txt in example subfolders with batch_CMakeLists.txt

for /R %%f in (CMakeLists.txt) do (
	if exist "%%f" (
		xcopy /qy "%~dp0batch_CMakeLists.txt" "%%f"
		echo replaced "%%f"
	) else (
	    rem file doesn't exist
	)
	
)

