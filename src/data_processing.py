"""
数据处理模块
处理序列数据的读取、存储和分析
"""

import os
import csv
import json
from typing import List, Dict, Any, Optional

class DataHandler:
    """
    数据处理器
    """
    
    @staticmethod
    def read_fasta(file_path: str) -> Dict[str, str]:
        """
        读取FASTA文件
        :return: {序列名: 序列}字典
        """
        sequences = {}
        current_seq = ""
        current_name = ""
        
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_name and current_seq:
                        sequences[current_name] = current_seq
                    current_name = line[1:]
                    current_seq = ""
                else:
                    current_seq += line
        
        if current_name and current_seq:
            sequences[current_name] = current_seq
        
        return sequences
    
    @staticmethod
    def write_fasta(sequences: Dict[str, str], file_path: str) -> None:
        """
        写入FASTA文件
        """
        with open(file_path, 'w') as f:
            for name, seq in sequences.items():
                f.write(f">{name}\n")
                # 每行60个字符
                for i in range(0, len(seq), 60):
                    f.write(seq[i:i+60] + "\n")
    
    @staticmethod
    def read_csv(file_path: str) -> List[Dict[str, str]]:
        """
        读取CSV文件
        :return: 字典列表
        """
        with open(file_path, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            return list(reader)
    
    @staticmethod
    def write_csv(data: List[Dict], file_path: str) -> None:
        """
        写入CSV文件
        """
        if not data:
            return
        
        with open(file_path, 'w', encoding='utf-8', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=data[0].keys())
            writer.writeheader()
            writer.writerows(data)
    
    @staticmethod
    def read_json(file_path: str) -> Any:
        """
        读取JSON文件
        """
        with open(file_path, 'r', encoding='utf-8') as f:
            return json.load(f)
    
    @staticmethod
    def write_json(data: Any, file_path: str) -> None:
        """
        写入JSON文件
        """
        with open(file_path, 'w', encoding='utf-8') as f:
            json.dump(data, f, ensure_ascii=False, indent=2)
    
    @staticmethod
    def load_target_genes(file_path: str) -> List[Dict]:
        """
        加载目标基因列表
        """
        return DataHandler.read_csv(file_path)
    
    @staticmethod
    def save_design_results(results: List[Dict], file_path: str) -> None:
        """
        保存设计结果
        """
        DataHandler.write_csv(results, file_path)


class SequenceDatabase:
    """
    序列数据库管理
    """
    
    def __init__(self, db_path: str = "data/sequences"):
        self.db_path = db_path
        self.use_memory = False
        
        try:
            os.makedirs(db_path, exist_ok=True)
        except (PermissionError, OSError):
            self.use_memory = True
            self.memory_db = {}
            
            import tempfile
            self.db_path = tempfile.mkdtemp(prefix="miRNA_designer_")
    
    def add_sequence(self, name: str, sequence: str, description: str = "") -> None:
        """
        添加序列到数据库
        """
        data = {
            'name': name,
            'sequence': sequence,
            'description': description,
            'length': len(sequence)
        }
        
        file_path = os.path.join(self.db_path, f"{name}.json")
        DataHandler.write_json(data, file_path)
    
    def get_sequence(self, name: str) -> Optional[Dict]:
        """
        获取序列
        """
        file_path = os.path.join(self.db_path, f"{name}.json")
        
        if os.path.exists(file_path):
            return DataHandler.read_json(file_path)
        
        return None
    
    def list_sequences(self) -> List[str]:
        """
        列出所有序列名称
        """
        sequences = []
        
        for filename in os.listdir(self.db_path):
            if filename.endswith('.json'):
                sequences.append(filename[:-5])
        
        return sequences
    
    def delete_sequence(self, name: str) -> bool:
        """
        删除序列
        """
        file_path = os.path.join(self.db_path, f"{name}.json")
        
        if os.path.exists(file_path):
            os.remove(file_path)
            return True
        
        return False