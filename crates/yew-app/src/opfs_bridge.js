async function getRoot() {
  if (!navigator.storage || !navigator.storage.getDirectory) {
    throw new Error("OPFS is not supported in this browser");
  }
  return navigator.storage.getDirectory();
}

async function getParentDir(root, parts) {
  let dir = root;
  for (const part of parts) {
    dir = await dir.getDirectoryHandle(part, { create: true });
  }
  return dir;
}

async function walk(dir, prefix = "") {
  const entries = [];
  for await (const [name, handle] of dir.entries()) {
    const fullPath = prefix ? `${prefix}/${name}` : name;
    if (handle.kind === "directory") {
      entries.push(...(await walk(handle, fullPath)));
    } else {
      const file = await handle.getFile();
      entries.push({
        path: fullPath,
        size: file.size,
        last_modified: file.lastModified,
      });
    }
  }
  return entries.sort((a, b) => a.path.localeCompare(b.path));
}

async function getFileFromPath(path) {
  const parts = path.split("/").filter(Boolean);
  const fileName = parts.pop();
  if (!fileName) {
    throw new Error("invalid OPFS path");
  }

  let dir = await getRoot();
  for (const part of parts) {
    dir = await dir.getDirectoryHandle(part);
  }

  const fileHandle = await dir.getFileHandle(fileName);
  return fileHandle.getFile();
}

globalThis.remOpfs = {
  async writeTextFile(path, content) {
    const parts = path.split("/").filter(Boolean);
    const fileName = parts.pop();
    if (!fileName) {
      throw new Error("invalid OPFS path");
    }

    const root = await getRoot();
    const parent = await getParentDir(root, parts);
    const handle = await parent.getFileHandle(fileName, { create: true });
    const writable = await handle.createWritable();
    await writable.write(content);
    await writable.close();
  },

  async listFiles() {
    const root = await getRoot();
    return walk(root);
  },

  async readTextFile(path) {
    const file = await getFileFromPath(path);
    return file.text();
  },

  async downloadTextFile(path) {
    const file = await getFileFromPath(path);
    const url = URL.createObjectURL(file);
    const anchor = document.createElement("a");
    anchor.href = url;
    anchor.download = path.split("/").pop() || "output.txt";
    document.body.appendChild(anchor);
    anchor.click();
    anchor.remove();
    URL.revokeObjectURL(url);
  },
};