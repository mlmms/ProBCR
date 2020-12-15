import os

src = "/Volumes/TOSHIBA_EXT/External_ProBCR"
count = 0
nrdirties = 0
strdirty = '_dirty'

parent_folders_dirty = []
new_parent_folders_dirty = []
child_folders_dirty = []
new_child_folders_dirty = []

for root, dirs, names in os.walk(src):
    #print(root)
    #print(dirs)
    #print(names)
    #count += 1
    #print("--------------------------------")
    count += 1
    if 'dirty' in names:
        nrdirties += 1

        if len(dirs) != 0:                  # is parent folder?
            parent_folders_dirty.append(root)
            new_parent_folders_dirty.append(root + strdirty)
            #print("parent_folders_dirty:", parent_folders_dirty)

        else:                               # is child folder?
            child_folders_dirty.append(root)
            new_child_folders_dirty.append(root + strdirty)

            # if child folder has "dirty", parent folder needs to be identified:
            parentstr = root[:root.rfind("/")]
            parent_folders_dirty.append(parentstr)
            new_parent_folders_dirty.append(parentstr + strdirty)

print("count = ", count)
print("nr dirties = ", nrdirties)
print("before parent_folders_dirty = ", parent_folders_dirty)
print("before child_folders_dirty = ", child_folders_dirty)

#remove duplicates: use 'set'
child_folders_dirty = list(set(child_folders_dirty))
parent_folders_dirty = list(set(parent_folders_dirty))
new_child_folders_dirty = list(set(new_child_folders_dirty))
new_parent_folders_dirty = list(set(new_parent_folders_dirty))

print("unique parent_folders_dirty = ", parent_folders_dirty)
print("unique child_folders_dirty = ", child_folders_dirty)

print("unique new parent_folders_dirty = ", new_parent_folders_dirty)
print("unique new child_folders_dirty = ", new_child_folders_dirty)


#renaming child folders (before renaming the "parent" folder!)
for child_folder, new_child_folder in zip(child_folders_dirty, new_child_folders_dirty):
    os.rename(child_folder, new_child_folder)

#renaming parent folders (before renaming the "parent" folder!)
for parent_folder, new_parent_folder in zip(parent_folders_dirty, new_parent_folders_dirty):
    os.rename(parent_folder, new_parent_folder)
