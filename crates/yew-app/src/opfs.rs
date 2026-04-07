use js_sys::Promise;
use serde::{Deserialize, Serialize};
use wasm_bindgen::prelude::*;
use wasm_bindgen_futures::JsFuture;

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct OpfsEntry {
    pub path: String,
    pub size: u64,
    pub last_modified: f64,
}

#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(js_namespace = remOpfs, js_name = writeTextFile)]
    fn write_text_file_js(path: &str, content: &str) -> Promise;

    #[wasm_bindgen(js_namespace = remOpfs, js_name = listFiles)]
    fn list_files_js() -> Promise;

    #[wasm_bindgen(js_namespace = remOpfs, js_name = readTextFile)]
    fn read_text_file_js(path: &str) -> Promise;

    #[wasm_bindgen(js_namespace = remOpfs, js_name = downloadTextFile)]
    fn download_text_file_js(path: &str) -> Promise;

    #[wasm_bindgen(js_namespace = remOpfs, js_name = deleteDir)]
    fn delete_dir_js(path: &str) -> Promise;
}

pub async fn write_text_file(path: &str, content: &str) -> Result<(), String> {
    JsFuture::from(write_text_file_js(path, content))
        .await
        .map(|_| ())
        .map_err(|err| format!("write_text_file failed: {:?}", err))
}

pub async fn list_files() -> Result<Vec<OpfsEntry>, String> {
    let value = JsFuture::from(list_files_js())
        .await
        .map_err(|err| format!("list_files failed: {:?}", err))?;
    serde_wasm_bindgen::from_value(value).map_err(|err| format!("decode list_files failed: {}", err))
}

pub async fn read_text_file(path: &str) -> Result<String, String> {
    let value = JsFuture::from(read_text_file_js(path))
        .await
        .map_err(|err| format!("read_text_file failed: {:?}", err))?;
    value.as_string().ok_or_else(|| "read_text_file returned non-string".to_string())
}

pub async fn download_text_file(path: &str) -> Result<(), String> {
    JsFuture::from(download_text_file_js(path))
        .await
        .map(|_| ())
        .map_err(|err| format!("download_text_file failed: {:?}", err))
}

pub async fn delete_dir(path: &str) -> Result<(), String> {
    JsFuture::from(delete_dir_js(path))
        .await
        .map(|_| ())
        .map_err(|err| format!("delete_dir failed: {:?}", err))
}