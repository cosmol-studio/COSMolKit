use serde::Deserialize;

#[derive(Deserialize)]
pub(super) struct Record {
    pub title: String,
    pub url: String,
    pub summary: String,
    text: String,
    #[serde(skip)]
    title_key: String,
}

pub(super) struct SearchIndex {
    pub records: Vec<Record>,
}

fn normalize(value: &str) -> String {
    value
        .chars()
        .flat_map(char::to_lowercase)
        .map(|c| if c.is_alphanumeric() { c } else { ' ' })
        .collect::<String>()
        .split_whitespace()
        .collect::<Vec<_>>()
        .join(" ")
}

impl SearchIndex {
    pub fn load(json: &str) -> Result<Self, serde_json::Error> {
        let mut records: Vec<Record> = serde_json::from_str(json)?;
        for record in &mut records {
            record.title_key = normalize(&record.title);
            record.text = normalize(&record.text);
        }
        Ok(Self { records })
    }

    pub fn search(&self, query: &str) -> Vec<usize> {
        let query = normalize(query);
        let terms: Vec<_> = query.split_whitespace().collect();
        if terms.is_empty() {
            return Vec::new();
        }
        let mut matches = Vec::new();
        for (id, record) in self.records.iter().enumerate() {
            if !terms
                .iter()
                .all(|term| record.title_key.contains(term) || record.text.contains(term))
            {
                continue;
            }
            let leaf = normalize(record.title.rsplit('.').next().unwrap_or(&record.title));
            let mut score = if leaf == query {
                300
            } else if record.title_key == query {
                250
            } else {
                0
            };
            if record.title_key.contains(&query) {
                score += 100;
            }
            score += terms
                .iter()
                .filter(|term| record.title_key.contains(**term))
                .count()
                * 20;
            matches.push((id, score));
        }
        matches.sort_by(|(a, sa), (b, sb)| {
            sb.cmp(sa)
                .then_with(|| self.records[*a].title.cmp(&self.records[*b].title))
                .then_with(|| self.records[*a].url.cmp(&self.records[*b].url))
        });
        matches.into_iter().map(|(id, _)| id).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn index() -> SearchIndex {
        SearchIndex::load(r##"[
            {"title":"Quick start","url":"/python/quickstart","summary":"Build molecules","text":"Molecule and atom examples, fingerprints"},
            {"title":"cosmolkit.Molecule.add_atom","url":"/python/api#cosmolkit.Molecule.add_atom","summary":"Python method","text":"Add an atom"},
            {"title":"cosmolkit.Molecule","url":"/python/api#cosmolkit.Molecule","summary":"Python class","text":"Molecular graph"},
            {"title":"分子结构","url":"/python/molecule","summary":"Unicode","text":"分子和原子"},
            {"title":"cosmolkit.from_smiles","url":"/python/api#from_smiles","summary":"Python function","text":"Create a molecule"}
        ]"##).unwrap()
    }

    #[test]
    fn exact_symbol_ranks_before_members_and_body_mentions() {
        let index = index();
        let results = index.search("Molecule");
        assert_eq!(results[0], 2);
        assert!(results.iter().position(|id| *id == 1) < results.iter().position(|id| *id == 0));
    }

    #[test]
    fn symbol_separators_case_and_partial_names_are_searchable() {
        let index = index();
        assert_eq!(index.search("FROM_SMILES"), vec![4]);
        assert_eq!(index.search("from smiles"), vec![4]);
        assert_eq!(index.search("fingerprint"), vec![0]);
        assert_eq!(index.records[4].url, "/python/api#from_smiles");
        assert_eq!(index.records[4].summary, "Python function");
    }

    #[test]
    fn every_query_term_is_required() {
        let index = index();
        assert_eq!(index.search("molecule & atom"), vec![1, 0]);
        assert!(index.search("molecule nonexistent").is_empty());
    }

    #[test]
    fn empty_punctuation_and_unknown_queries_have_no_results() {
        for query in ["", " \t", "---", "<script>alert(1)</script>"] {
            assert!(index().search(query).is_empty());
        }
    }

    #[test]
    fn unicode_search_and_order_are_stable() {
        let index = index();
        assert_eq!(index.search("分子"), vec![3]);
        assert_eq!(index.search("Molecule"), index.search("molecule"));
    }

    #[test]
    fn malformed_index_is_an_error() {
        assert!(SearchIndex::load("not JSON").is_err());
        assert!(SearchIndex::load(r#"[{"title":"missing fields"}]"#).is_err());
    }
}
