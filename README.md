# uHAF: Universal Hierarchical Annotation Framework

uHAF is a Python library developed to address the challenges of inconsistent cell type annotations in single-cell transcriptomics, such as varied naming conventions and hierarchical granularity. It integrates organ-specific hierarchical cell type trees (uHAF-T) and a mapping tool (uHAF-GPT) powered by large language models to provide a standardized framework for annotation. By enabling consistent label unification, hierarchical analysis, and integration of diverse datasets, uHAF enhances machine learning applications and facilitates biologically meaningful evaluations. This library is an essential resource for the single-cell research community, fostering collaborative refinement and standardization of cell type annotations.

## Explore Online

- [uHAF-T Explorer](https://uhaf.unifiedcellatlas.org): Browse and explore uHAF-Ts.
- [uHAF Mapping](https://uhaf.unifiedcellatlas.org/#/uHAFMapping): Map custom cell type labels to uHAF-Ts.

## Installation

Install uHAF via pip:

```bash
pip install uhaf
```

## Getting Started

### Building uHAF

Start by building a uHAF object for your dataset:

```python
import uhaf as uhaflib

uhaf = uhaflib.build_uhaf('uHAF2.2.0')
print(len(uhaf.df_uhafs))
```

This generates a uHAF instance containing annotations for all organs. The example above initializes the `uHAF2.2.0` dataset.

### Tracing Cell Types

Trace the hierarchical ancestry of a target cell type:

```python
ancestors = uhaf.track_cell_from_uHAF(sheet_name='Lung', cell_type_target='CD8 T cell')
print(ancestors)
```

Output:

```
['Cell', 'Lymphocyte', 'T cell', 'CD8 T cell']
```

### Annotation Levels

Retrieve hierarchical annotation levels for cell types. Specify the desired level (e.g., main, middle, or fine).

```python
example_cell_types = ['Pericyte', 'Macrophage', 'Monocyte-derived macrophage', 'Monocyte', 'Dendritic cell']
annotation_level = 2  # Middle cell type level
annotations = uhaf.set_annotation_level(example_cell_types, sheet_name='Heart', annotation_level=annotation_level)
print(annotations)
```

Example Output:

```
{'Pericyte': 'Pericyte', 'Macrophage': 'Macrophage', 'Monocyte-derived macrophage': 'Macrophage', 'Monocyte': 'Monocyte', 'Dendritic cell': 'Dendritic cell'}
```

### Mapping Custom Labels

To map custom cell type labels to uHAF:

1. Prepare unique cell type labels from your dataset:

   ```python
   original_labels = ['V-CM', 'LA-CM', 'RA-CM', 'Capillary-EC', 'Lymphatic-EC']
   ```

2. Generate uHAF-GPT prompts:

   ```python
   print(uhaf.generate_uhaf_GPTs_prompts('Heart', original_labels))
   ```

   Copy the output and use it on the [uHAF-GPT Mapping Website](https://uhaf.unifiedcellatlas.org/#/uHAFMapping) to get the mapped labels.

3. Use the mapping dictionary to transform your labels:

   ```python
   mapping_results = {"V-CM": "Ventricle cardiomyocyte cell", "LA-CM": "Atrial cardiomyocyte cell"}
   transformed_labels = [mapping_results[label] for label in original_labels]
   print(transformed_labels)
   ```

### Generating Nested JSON

Export the uHAF tree for a specific organ in nested JSON format:

```python
print(uhaf.dict_uhafs['Heart'])
```

## Contribution

We welcome contributions to improve and expand the uHAF library. For more details, please refer to our contribution guidelines.

## License

This project is licensed under the MIT License. See the LICENSE file for details.
