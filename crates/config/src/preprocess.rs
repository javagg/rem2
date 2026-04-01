/// Strip C++-style `//` line comments and `/* */` block comments from JSON-like text,
/// while leaving string literals untouched.
pub fn strip_comments(s: &str) -> String {
    let bytes = s.as_bytes();
    let n = bytes.len();
    let mut out = String::with_capacity(n);
    let mut i = 0;
    let mut in_string = false;

    while i < n {
        let b = bytes[i];

        if in_string {
            out.push(b as char);
            if b == b'\\' && i + 1 < n {
                // escaped character — emit both and skip ahead
                i += 1;
                out.push(bytes[i] as char);
            } else if b == b'"' {
                in_string = false;
            }
            i += 1;
            continue;
        }

        // check for comment start
        if b == b'/' && i + 1 < n {
            let next = bytes[i + 1];
            if next == b'/' {
                // line comment — skip until newline
                while i < n && bytes[i] != b'\n' {
                    i += 1;
                }
                continue;
            }
            if next == b'*' {
                // block comment — skip until */
                i += 2;
                while i + 1 < n {
                    if bytes[i] == b'*' && bytes[i + 1] == b'/' {
                        i += 2;
                        break;
                    }
                    i += 1;
                }
                continue;
            }
        }

        if b == b'"' {
            in_string = true;
        }
        out.push(b as char);
        i += 1;
    }
    out
}

/// Parse a Palace attribute spec into a sorted Vec<u32>.
/// Accepts either a JSON array `[1,2,3]` (already deserialized) OR
/// a string like `"1,3-5,7"` → `[1,3,4,5,7]`.
pub fn expand_ranges(s: &str) -> Result<Vec<u32>, String> {
    let s = s.trim();
    let mut result = Vec::new();

    for part in s.split(',') {
        let part = part.trim();
        if part.is_empty() {
            continue;
        }
        if let Some(dash) = part.find('-') {
            let lo: u32 = part[..dash]
                .trim()
                .parse()
                .map_err(|_| format!("invalid range start: '{}'", &part[..dash]))?;
            let hi: u32 = part[dash + 1..]
                .trim()
                .parse()
                .map_err(|_| format!("invalid range end: '{}'", &part[dash + 1..]))?;
            if lo > hi {
                return Err(format!("range {}-{} is inverted", lo, hi));
            }
            result.extend(lo..=hi);
        } else {
            let v: u32 = part
                .parse()
                .map_err(|_| format!("invalid integer: '{}'", part))?;
            result.push(v);
        }
    }
    result.sort_unstable();
    result.dedup();
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn strip_line_comment() {
        let s = r#"{ "a": 1 // comment
}"#;
        let out = strip_comments(s);
        assert!(out.contains("\"a\": 1"));
        assert!(!out.contains("comment"));
    }

    #[test]
    fn strip_block_comment() {
        let s = r#"{ /* block */ "b": 2 }"#;
        let out = strip_comments(s);
        assert!(out.contains("\"b\": 2"));
        assert!(!out.contains("block"));
    }

    #[test]
    fn preserve_string_with_slashes() {
        let s = r#"{ "url": "http://example.com" }"#;
        let out = strip_comments(s);
        assert!(out.contains("http://example.com"));
    }

    #[test]
    fn expand_simple() {
        assert_eq!(expand_ranges("1,2,3").unwrap(), vec![1, 2, 3]);
    }

    #[test]
    fn expand_range() {
        assert_eq!(expand_ranges("1,3-5,7").unwrap(), vec![1, 3, 4, 5, 7]);
    }

    #[test]
    fn expand_dedup() {
        assert_eq!(expand_ranges("1,1,2").unwrap(), vec![1, 2]);
    }
}
