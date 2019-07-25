;
; Routines for reading Peano-Hilbert key sorted Eagle snapshots
;


; Convenience function for reading full datasets
function re_read_dataset, file_id, name
dset_id = h5d_open(file_id, name)
data    = h5d_read(dset_id)
h5d_close, dset_id
return, data
end

; Convenience function for reading attributes of groups
function re_read_attribute, file_id, group_name, attr_name
group_id = H5G_open(file_id, group_name)
attr_id  = H5A_open_name(group_id, attr_name)
data = h5a_read(attr_id)
h5a_close, attr_id
h5g_close, group_id
return, data
end

function str, s
return, strtrim(string(s),2)
end

;
; Function to calculate P-H keys
;
; To get adequate performance out of this we need to calculate
; a large array of keys all in one go.
;
function peano_hilbert_keys, ix, iy, iz, bits

rottable3 = [ $
  [36, 28, 25, 27, 10, 10, 25, 27], $
  [29, 11, 24, 24, 37, 11, 26, 26], $
  [8, 8, 25, 27, 30, 38, 25, 27], $
  [9, 39, 24, 24, 9, 31, 26, 26], $
  [40, 24, 44, 32, 40, 6, 44, 6], $
  [25, 7, 33, 7, 41, 41, 45, 45], $
  [4, 42, 4, 46, 26, 42, 34, 46], $
  [43, 43, 47, 47, 5, 27, 5, 35], $
  [33, 35, 36, 28, 33, 35, 2, 2], $
  [32, 32, 29, 3, 34, 34, 37, 3], $
  [33, 35, 0, 0, 33, 35, 30, 38], $
  [32, 32, 1, 39, 34, 34, 1, 31], $
  [24, 42, 32, 46, 14, 42, 14, 46], $
  [43, 43, 47, 47, 25, 15, 33, 15], $
  [40, 12, 44, 12, 40, 26, 44, 34], $
  [13, 27, 13, 35, 41, 41, 45, 45], $
  [28, 41, 28, 22, 38, 43, 38, 22], $
  [42, 40, 23, 23, 29, 39, 29, 39], $
  [41, 36, 20, 36, 43, 30, 20, 30], $
  [37, 31, 37, 31, 42, 40, 21, 21], $
  [28, 18, 28, 45, 38, 18, 38, 47], $
  [19, 19, 46, 44, 29, 39, 29, 39], $
  [16, 36, 45, 36, 16, 30, 47, 30], $
  [37, 31, 37, 31, 17, 17, 46, 44], $
  [12, 4, 1, 3, 34, 34, 1, 3], $
  [5, 35, 0, 0, 13, 35, 2, 2], $
  [32, 32, 1, 3, 6, 14, 1, 3], $
  [33, 15, 0, 0, 33, 7, 2, 2], $
  [16, 0, 20, 8, 16, 30, 20, 30], $
  [1, 31, 9, 31, 17, 17, 21, 21], $
  [28, 18, 28, 22, 2, 18, 10, 22], $
  [19, 19, 23, 23, 29, 3, 29, 11], $
  [9, 11, 12, 4, 9, 11, 26, 26], $
  [8, 8, 5, 27, 10, 10, 13, 27], $
  [9, 11, 24, 24, 9, 11, 6, 14], $
  [8, 8, 25, 15, 10, 10, 25, 7], $
  [0, 18, 8, 22, 38, 18, 38, 22], $
  [19, 19, 23, 23, 1, 39, 9, 39], $
  [16, 36, 20, 36, 16, 2, 20, 10], $
  [37, 3, 37, 11, 17, 17, 21, 21], $
  [4, 17, 4, 46, 14, 19, 14, 46], $
  [18, 16, 47, 47, 5, 15, 5, 15], $
  [17, 12, 44, 12, 19, 6, 44, 6], $
  [13, 7, 13, 7, 18, 16, 45, 45], $
  [4, 42, 4, 21, 14, 42, 14, 23], $
  [43, 43, 22, 20, 5, 15, 5, 15], $
  [40, 12, 21, 12, 40, 6, 23, 6], $
  [13, 7, 13, 7, 41, 41, 22, 20]]

subpix3 = [ $
  [0, 7, 1, 6, 3, 4, 2, 5], $
  [7, 4, 6, 5, 0, 3, 1, 2], $
  [4, 3, 5, 2, 7, 0, 6, 1], $
  [3, 0, 2, 1, 4, 7, 5, 6], $
  [1, 0, 6, 7, 2, 3, 5, 4], $
  [0, 3, 7, 4, 1, 2, 6, 5], $
  [3, 2, 4, 5, 0, 1, 7, 6], $
  [2, 1, 5, 6, 3, 0, 4, 7], $
  [6, 1, 7, 0, 5, 2, 4, 3], $
  [1, 2, 0, 3, 6, 5, 7, 4], $
  [2, 5, 3, 4, 1, 6, 0, 7], $
  [5, 6, 4, 7, 2, 1, 3, 0], $
  [7, 6, 0, 1, 4, 5, 3, 2], $
  [6, 5, 1, 2, 7, 4, 0, 3], $
  [5, 4, 2, 3, 6, 7, 1, 0], $
  [4, 7, 3, 0, 5, 6, 2, 1], $
  [6, 7, 5, 4, 1, 0, 2, 3], $
  [7, 0, 4, 3, 6, 1, 5, 2], $
  [0, 1, 3, 2, 7, 6, 4, 5], $
  [1, 6, 2, 5, 0, 7, 3, 4], $
  [2, 3, 1, 0, 5, 4, 6, 7], $
  [3, 4, 0, 7, 2, 5, 1, 6], $
  [4, 5, 7, 6, 3, 2, 0, 1], $
  [5, 2, 6, 1, 4, 3, 7, 0], $
  [7, 0, 6, 1, 4, 3, 5, 2], $
  [0, 3, 1, 2, 7, 4, 6, 5], $
  [3, 4, 2, 5, 0, 7, 1, 6], $
  [4, 7, 5, 6, 3, 0, 2, 1], $
  [6, 7, 1, 0, 5, 4, 2, 3], $
  [7, 4, 0, 3, 6, 5, 1, 2], $
  [4, 5, 3, 2, 7, 6, 0, 1], $
  [5, 6, 2, 1, 4, 7, 3, 0], $
  [1, 6, 0, 7, 2, 5, 3, 4], $
  [6, 5, 7, 4, 1, 2, 0, 3], $
  [5, 2, 4, 3, 6, 1, 7, 0], $
  [2, 1, 3, 0, 5, 6, 4, 7], $
  [0, 1, 7, 6, 3, 2, 4, 5], $
  [1, 2, 6, 5, 0, 3, 7, 4], $
  [2, 3, 5, 4, 1, 0, 6, 7], $
  [3, 0, 4, 7, 2, 1, 5, 6], $
  [1, 0, 2, 3, 6, 7, 5, 4], $
  [0, 7, 3, 4, 1, 6, 2, 5], $
  [7, 6, 4, 5, 0, 1, 3, 2], $
  [6, 1, 5, 2, 7, 0, 4, 3], $
  [5, 4, 6, 7, 2, 3, 1, 0], $
  [4, 3, 7, 0, 5, 2, 6, 1], $
  [3, 2, 0, 1, 4, 5, 7, 6], $
  [2, 5, 1, 6, 3, 4, 0, 7]]

n        = long(n_elements(ix))
key      = lon64arr(n)
rotation = 0L

for i = bits-1, 0, -1 do begin
    mask = ishft(1L, i)
    
    px = intarr(n)
    ind = where(ix and mask)
    if ind[0] ne -1 then px[ind] = 4

    py = intarr(n)
    ind = where(iy and mask)
    if ind[0] ne -1 then py[ind] = 2

    pz = intarr(n)
    ind = where(iz and mask)
    if ind[0] ne -1 then pz[ind] = 1

    pix = px or py or pz

    key = ishft(key, 3)
    key = key or subpix3[pix, rotation]
    rotation = rottable3[pix, rotation]

endfor

return, key
end

;
; Open a snapshot by specifying the name of one file from it
;
; Returns a struct with all the necessary information to read
; sub-regions of the snapshot.
;
function open_snapshot, fname

; Read header stuff we need
file_id  = h5f_open(fname)
boxsize  = re_read_attribute(file_id, "Header", "BoxSize")
numfiles = re_read_attribute(file_id, "Header", "NumFilesPerSnapshot")
nptot    = re_read_attribute(file_id, "Header", "NumPart_Total")
nptot_hw = re_read_attribute(file_id, "Header", "NumPart_Total_HighWord")
nptot    = long64(nptot) + long64(nptot_hw) * 2LL^32
hashbits = re_read_attribute(file_id, "HashTable", "HashBits")

; Find size of hash map, allocate and initialise
ncell    = 2L^hashbits
nhash    = ncell^3
hashmap  = intarr(nhash)
hashmap[*] = 0

; Find snapshot base name
basename = stregex(fname, "(.*)\.[0-9]+\.hdf5", /extract, /subexpr)
basename = basename[1]

; Allocate and read in range of keys in each file
first_key_in_file = lon64arr(6, numfiles)
last_key_in_file  = lon64arr(6, numfiles)
num_keys_in_file  = lon64arr(6, numfiles)
for i = 0, 5, 1 do begin
    if nptot[i] gt 0 then begin
        first_key_in_file[i,*] = re_read_dataset(file_id, $
                                                 "HashTable/PartType"+ $
                                                 str(i)+"/FirstKeyInFile")
        last_key_in_file[i,*] = re_read_dataset(file_id, $
                                                "HashTable/PartType"+ $
                                                str(i)+"/LastKeyInFile")
        num_keys_in_file[i,*] = re_read_dataset(file_id, $
                                                "HashTable/PartType"+ $
                                                str(i)+"/NumKeysInFile")
    endif
endfor

; Set up pointers to hash table data which we'll read as needed later
part_per_cell = ptrarr(6, numfiles)
first_in_cell = ptrarr(6, numfiles)

; Determine which file each hash cell is in
filemap    = intarr(6,nhash)
filemap[*,*] = -1 ; Certain (empty) cells may be in no file 
for itype=0, 5, 1 do begin
    if nptot[itype] gt 0 then begin
        for ifile=0, numfiles-1, 1 do begin
            filemap[itype, first_key_in_file[itype, ifile]: $
                    last_key_in_file[itype, ifile]] = ifile
        endfor
    endif
endfor

snap = create_struct( $
                      "boxsize",           boxsize, $
                      "numfiles",          numfiles, $
                      "nptot",             nptot, $
                      "hashbits",          hashbits, $
                      "ncell",             ncell, $
                      "nhash",             nhash, $
                      "hashmap",           hashmap, $
                      "basename",          basename, $
                      "first_key_in_file", first_key_in_file, $
                      "last_key_in_file",  last_key_in_file, $
                      "num_keys_in_file",  num_keys_in_file, $
                      "part_per_cell",     part_per_cell, $
                      "first_in_cell",     first_in_cell, $
                      "filemap",           filemap $
                    )

return, snap
end


;
; Close a snapshot
;
; Deallocates any hash table information that was read in
;
pro close_snapshot, snap

for i = 0, 5, 1 do begin
    for j = 0, snap.numfiles-1, 1 do begin
        if (ptr_valid(snap.part_per_cell[i,j])) then begin
            ptr_free, snap.part_per_cell[i,j]
        endif
    endfor
endfor

end


;
; Select a region to read
;
pro select_region, snap, xmin, xmax, ymin, ymax, zmin, zmax

ixmin = floor(xmin / snap.boxsize * snap.ncell)
ixmax = floor(xmax / snap.boxsize * snap.ncell)
iymin = floor(ymin / snap.boxsize * snap.ncell)
iymax = floor(ymax / snap.boxsize * snap.ncell)
izmin = floor(zmin / snap.boxsize * snap.ncell)
izmax = floor(zmax / snap.boxsize * snap.ncell)

nx = ixmax - ixmin + 1
ny = iymax - iymin + 1
nz = izmax - izmin + 1

coords = array_indices([nx,ny,nz], indgen(nx*ny*nz), /dimensions)

if nx*ny*nz gt 1 then begin
    ix = coords[0,*] + ixmin
    iy = coords[1,*] + iymin
    iz = coords[2,*] + izmin
endif else begin
    ix = [coords[0] + ixmin]
    iy = [coords[1] + iymin]
    iz = [coords[2] + izmin]
endelse

keys = peano_hilbert_keys(ix,iy,iz,snap.hashbits)

snap.hashmap[keys] = 1

end


;
; Clear currently selected region
;
pro clear_selection, snap

snap.hashmap[*] = 0

end

;
; Load the hash table for one type in one file
;
pro load_hash_table, snap, itype, ifile

; Check if already loaded
if ptr_valid(snap.part_per_cell[itype, ifile]) then return

; Read in the dataset
dset_name = "HashTable/PartType"+str(itype)+"/NumParticleInCell"
fname = snap.basename + "." + str(ifile) + ".hdf5"
file_id = h5f_open(fname)
data = re_read_dataset(file_id, dset_name)
h5f_close, file_id
snap.part_per_cell[itype, ifile] = ptr_new(data, /no_copy)

; Calculate offset to start of each cell
data = total(*snap.part_per_cell[itype, ifile], /cumulative, /integer) - $
  *snap.part_per_cell[itype, ifile]

snap.first_in_cell[itype, ifile] = ptr_new(data, /no_copy)

end

;
; Count selected particles
;
function count_particles, snap, itype

; Check if there are any particles of this type
if snap.nptot[itype] eq 0 then begin
    return, 0
endif

np = 0LL
for ifile=0, snap.numfiles-1, 1 do begin

    if snap.num_keys_in_file[itype, ifile] gt 0 then begin

        ; Check if we need to read from this file
        local_hashmap = snap.hashmap[snap.first_key_in_file[itype,ifile]: $
                                     snap.last_key_in_file[itype,ifile]]
        n = total(local_hashmap, /integer)

        if n gt 0 then begin
            ; Need to read from this one. Make sure we have the hash table.
            load_hash_table, snap, itype, ifile
            ; Count particles to read in this file
            ppc = *(snap.part_per_cell[itype, ifile])
            ind = where((local_hashmap ne 0) and (ppc gt 0), ncell_read)
            if ncell_read gt 0 then np = np + total(ppc[ind], /integer)
        endif

    endif

endfor

return, np
end

;
; Determine type and rank of a dataset
;
pro get_dataset_info, snap, itype, dset_name, typecode, rank

; Check we have some of these particles
if snap.nptot[itype] eq 0 then begin
    print, "There are no particles of the specified type!"
    stop
endif

; Find a file with particles of this type
ifile = 0
while snap.num_keys_in_file[itype, ifile] eq 0 do begin
    ifile = ifile + 1
endwhile

; Open this file
fname = snap.basename + "." + str(ifile) + ".hdf5"
file_id = h5f_open(fname)

; Open the dataset
name = "PartType"+str(itype)+"/"+dset_name
dset_id   = h5d_open(file_id, name)
dspace_id = h5d_get_space(dset_id)
dtype_id  = h5d_get_type(dset_id)
rank  = h5s_get_simple_extent_ndims(dspace_id)
class = h5t_get_class(dtype_id) 
size  = h5t_get_size(dtype_id)
h5t_close, dtype_id
h5s_close, dspace_id
h5d_close, dset_id
h5f_close, file_id

typecode = -1
if class eq "H5T_INTEGER" then begin
    if size le 4 then begin
        typecode = 0
    endif else begin
        typecode = 1
    endelse
endif
if class eq "H5T_FLOAT" then begin
    if size le 4 then begin
        typecode = 2
    endif else begin
        typecode = 3
    endelse
endif

end

;
; Read a dataset for the selected particles
;
function read_dataset, snap, itype, dset_name

; Check we have some of these particles
if snap.nptot[itype] eq 0 then begin
    print, "There are no particles of the specified type!"
    stop
endif

; Get size and type of result
n = count_particles(snap, itype)

; Check there are particles
if n eq 0 then begin
    print, "There are no particles of the specified type in the selected region!"
    stop
endif

get_dataset_info, snap, itype, dset_name, typecode, rank

; Set up output array and corresponding dataspace
if rank eq 1 then begin
    if typecode eq 0 then data = lonarr(n)
    if typecode eq 1 then data = lon64arr(n)
    if typecode eq 2 then data = fltarr(n)
    if typecode eq 3 then data = dblarr(n)
endif else begin
    if typecode eq 0 then data = lonarr(3,n)
    if typecode eq 1 then data = lon64arr(3,n)
    if typecode eq 2 then data = fltarr(3,n)
    if typecode eq 3 then data = dblarr(3,n)
endelse

np = 0LL

; Loop over all files
for ifile=0, snap.numfiles-1, 1 do begin
    
    if snap.num_keys_in_file[itype, ifile] gt 0 then begin

                                ; Find which hash cells we need from this file 
        local_hashmap = snap.hashmap[snap.first_key_in_file[itype,ifile]: $
                                     snap.last_key_in_file[itype,ifile]]

                                ; Check if we need to read from this file
        n = total(local_hashmap, /integer)

        if n gt 0 then begin

                                ; Need to read from this one. Make sure we have the hash table.
            load_hash_table, snap, itype, ifile

                                ; Count particles to read in this file
            ppc = *(snap.part_per_cell[itype, ifile])
            fic = *(snap.first_in_cell[itype, ifile])
            ind = where((local_hashmap ne 0) and (ppc gt 0), ncell_read)
            if ncell_read eq 0 then continue

            np_file = total(ppc[ind], /integer)

                                ; Find offsets and lengths of sections to read
            offsets = fic[ind]
            lengths = ppc[ind]

                                ; Open the file
            name = snap.basename+"."+str(ifile)+".hdf5"
            file_id = h5f_open(name)

                                ; Open the dataset and get its dataspace
            name = "PartType"+str(itype)+"/"+dset_name
            dset_id = h5d_open(file_id, name)
            dspace_id = h5d_get_space(dset_id)
            h5s_select_none, dspace_id

                                ; Loop over sections to read
            start = lonarr(2)
            count = lonarr(2)
            count[0] = 3 ; In case dataset is 2D, will be overwritten if 1D
            if rank eq 1 then begin
                                ; Scalar quantity
                n = 0LL
                for i=0LL, n_elements(lengths)-1, 1 do begin
                    start[0] = offsets[i]
                    count[0] = lengths[i]
                    n = n + count[0]
                    h5s_select_hyperslab, dspace_id, start[0:0], count[0:0]
                endfor
            endif else begin
                                ; Vector quantity
                for i=0LL, n_elements(lengths)-1, 1 do begin
                    start[1] = offsets[i]
                    count[1] = lengths[i]
                    h5s_select_hyperslab, dspace_id, start, count
                endfor
            endelse

                                ; Create memory dataspace
            if rank eq 1 then begin
                memspace_id = h5s_create_simple([np_file])
            endif else begin
                memspace_id = h5s_create_simple([3, np_file])
            endelse

                                ; Read the data
            if rank eq 1 then begin
                data[np:np+np_file-1] = h5d_read(dset_id, $
                                                 file_space=dspace_id, $
                                                 memory_space=memspace_id)
            endif else begin
                data[0:2,np:np+np_file-1] = h5d_read(dset_id, $
                                                     file_space=dspace_id, $
                                                     memory_space=memspace_id)  
            endelse

                                ; Close dataset and file etc
            h5s_close, dspace_id
            h5s_close, memspace_id
            h5d_close, dset_id
            h5f_close, file_id

            np = np + np_file
        endif

    endif

endfor

return, data
end






