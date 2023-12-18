import argparse
import json
import os
import re
import sys
from typing import Any, Dict
import WDL
from ruamel.yaml import YAML


def from_camel_case_to_snake_case(name: str) -> str:
    """Converts a camelCase string to snake_case."""
    name = re.sub('(.)([A-Z][a-z]+)', r'\1_\2', name)
    return re.sub('([a-z0-9])([A-Z])', r'\1_\2', name).lower()


def generate_inputs(inputs: Dict[str, Any], test_path: str) -> str:
    """Creates a JSON file with the provided inputs in the specified directory."""
    save_path = f'{test_path}/inputs.json'
    with open(save_path, 'w') as outfile:
        json.dump(inputs, outfile, indent=4)
    return save_path


def generate_yaml(task_name: str, tags: list[str], inputs_path: str, wdl_path: str, output_paths: list[str], test_path: str) -> None:
    """Creates a YAML file with test information for the provided task."""
    data = [{
        'name': f'Check if {task_name} produces the expected outputs',
        'tags': tags,
        'command': f'miniwdl run --task {task_name} -i {inputs_path} {wdl_path}',
        'files': [{'path': path} for path in output_paths]
    }]

    yaml = YAML()
    with open(f'{test_path}/test_{from_camel_case_to_snake_case(task_name)}_task.yaml', 'w') as outfile:
        yaml.dump(data, outfile)


def generate_tests_for_task(task: WDL.Task, task_path: str, test_base_path: str) -> None:
    """Generates test inputs and YAML files for a provided WDL task."""
    test_path = f"tests/{os.path.splitext(task_path)[0]}"

    if not os.path.exists(test_path):
        os.makedirs(test_path)

    inputs_path = generate_inputs({j.name: str(j.value.type) for j in task.required_inputs}, test_path)
    output_paths = [f'_LAST/out/{o.name}/<{o.info}>' for o in task.effective_outputs]

    tag = "/".join(test_path.split("/")[-2:])  # "tests/tasks/minimap2/align" -> "minimap2/align"
    generate_yaml(
        task_name=task.name,
        tags=[tag],
        inputs_path=inputs_path,
        wdl_path=task_path,
        output_paths=output_paths,
        test_path=test_path,
    )


def main() -> None:
    """Main script function. Processes WDL file and generates tests for the specified task."""
    parser = argparse.ArgumentParser(description="Generate test files for a specific WDL task.")
    parser.add_argument("--wdl_file_path", help="The path to the WDL file.")
    parser.add_argument("--task_name", help="The name of the task in the WDL file for which to generate test files.")
    parser.add_argument("--test_path", default="tasks/tests", help="The path to the test folder.")
    args = parser.parse_args()

    doc = WDL.load(args.wdl_file_path)
    for i in doc.tasks:
        if isinstance(i, WDL.Task) and i.name == args.task_name:
            print(f'Generating test for {i.name}')
            generate_tests_for_task(i, task_path=args.wdl_file_path, test_base_path=args.test_path)


if __name__ == "__main__":
    main()
