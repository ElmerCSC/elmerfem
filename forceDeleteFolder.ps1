$folderPath = ".\build"

try {
    Remove-Item -Path $folderPath -Recurse -Force -ErrorAction Stop
}
catch {
    # Assuming we have a custom function 'Get-LockingProcess' that can find the process locking a file
    $lockingProcess = Get-LockingProcess -Path $folderPath

    if ($lockingProcess) {
        Stop-Process -Id $lockingProcess.Id -Force
        Remove-Item -Path $folderPath -Recurse -Force
    }
    else {
        Write-Error "Unable to find the locking process."
    }
}
