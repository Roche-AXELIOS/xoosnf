def germlineLibraryTypeOption(runType) {
    def normalizedRunType = runType?.toString()?.toLowerCase()
    def libraryType = normalizedRunType in ['ys', 'ysu', 'ys-new'] ? 'sbx-s' : normalizedRunType
    if (!(libraryType in ['sbx-d', 'sbx-fast', 'sbx-s'])) {
        throw new IllegalArgumentException("Unsupported SVC germline run type '${runType}'. Expected SBX-D, SBX-FAST, YS, YSU, or YS-NEW.")
    }
    return "--library-type ${libraryType}"
}
