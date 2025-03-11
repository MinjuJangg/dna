import os
import pandas as pd
import argparse

def blastp_filtering(blastp_dir, output_list, keyword):
    """
    Filters BLASTP results to identify valid FASTA files based on a keyword and saves their names to a .txt file.

    Args:
        blastp_dir (str): Directory containing BLASTP result files.
        output_list (str): Path to save the list of valid FASTA file names.
        keyword (str): Keyword to filter results in 'salltitles' column.

    Returns:
        None
    """
    valid_fasta_files = set()  # 유효한 FASTA 파일명을 저장할 집합

    # BLASTP 결과 디렉토리 내 모든 파일 확인
    for filename in sorted(os.listdir(blastp_dir)):
        if filename.endswith("_result.txt"):  # BLASTP 결과 파일 처리
            input_file_path = os.path.join(blastp_dir, filename)
            try:
                # BLASTP 결과 읽기
                df = pd.read_csv(input_file_path, sep="\t")

                # 'salltitles' 열에서 지정된 keyword만 포함된 행이 있는지 필터링
                if df['salltitles'].apply(lambda x: all(keyword.lower() in item.lower() for item in x.split("<>"))).all():
                    # 유효한 FASTA 파일 이름 추출
                    fasta_filename = filename.replace("_result.txt", ".fa")
                    valid_fasta_files.add(fasta_filename)
                    print(f"Valid BLASTP result identified: {filename}")
                else:
                    print(f"File skipped (does not match keyword '{keyword}'): {filename}")
            except Exception as e:
                print(f"Error processing {filename}: {e}")

    # 유효한 FASTA 파일 이름 목록을 .txt 파일로 저장
    with open(output_list, "w") as file:
        for fasta_file in sorted(valid_fasta_files):
            file.write(fasta_file + "\n")
            print(f"Added to list: {fasta_file}")

    print(f"Valid FASTA file list saved to: {output_list}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter BLASTP results and save valid FASTA file names.")
    parser.add_argument("--blastp_dir", type=str, required=True, help="Directory containing BLASTP result files.")
    parser.add_argument("--output_list", type=str, required=True, help="Path to save the list of valid FASTA file names.")
    parser.add_argument("--keyword", type=str, required=True, help="Keyword to filter results in 'salltitles' column.")

    args = parser.parse_args()

    blastp_filtering(args.blastp_dir, args.output_list, args.keyword)
