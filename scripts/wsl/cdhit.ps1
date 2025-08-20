param(
    [Parameter()]$which = '',
    [Parameter(ValueFromRemainingArguments)]$cdhitArgs
)

process {
    echo $which $cdhitArgs
    $posixArgs = $cdhitArgs -split ' ' | ForEach-Object {
        if ($_ -match '\\') {
            $abs = (Resolve-Path $_).Path
            $fixes = (cygpath $abs).Trim('".')
            "/mnt$fixes"
        } else { $_ }
    }
    $cd_exe = "cd-hit"
    if ($which) { $cd_exe = "$cd_exe-$which"}
    echo $cd_exe @posixArgs
}