import requests
import time
import xml.etree.ElementTree as ET
from typing import Dict, List, Optional, Tuple
import json

class NCBIBLAST:
    """NCBI BLAST API客户端"""
    
    def __init__(self):
        self.base_url = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
    
    def submit_blast(self, sequence: str, program: str = "blastn", database: str = "refseq_rna", 
                     word_size: int = 7, expect: float = 10.0, hitlist_size: int = 50) -> Optional[str]:
        """
        提交BLAST请求
        
        参数:
            sequence: 查询序列
            program: BLAST程序类型 (blastn, megablast等)
            database: 目标数据库
            word_size: 种子词大小
            expect: E-value阈值
            hitlist_size: 返回结果数量
        
        返回:
            request_id: BLAST任务ID，如果失败返回None
        """
        try:
            params = {
                'CMD': 'Put',
                'PROGRAM': program,
                'DATABASE': database,
                'QUERY': sequence,
                'WORD_SIZE': word_size,
                'EXPECT': expect,
                'HITLIST_SIZE': hitlist_size,
                'FORMAT_TYPE': 'XML'
            }
            
            response = requests.post(self.base_url, data=params, timeout=30)
            
            if response.status_code == 200:
                # 从响应中提取request_id
                text = response.text
                start = text.find('RID = ')
                if start != -1:
                    end = text.find('\n', start)
                    if end != -1:
                        return text[start+6:end].strip()
            
            return None
        except Exception as e:
            print(f"BLAST提交失败: {e}")
            return None
    
    def check_status(self, request_id: str) -> Tuple[str, Optional[str]]:
        """
        检查BLAST任务状态
        
        返回:
            status: 任务状态 (WAITING, RUNNING, COMPLETE, FAILED)
            result: 结果XML字符串（如果完成）
        """
        try:
            params = {
                'CMD': 'Get',
                'RID': request_id,
                'FORMAT_TYPE': 'XML'
            }
            
            response = requests.get(self.base_url, params=params, timeout=30)
            
            if response.status_code == 200:
                text = response.text
                
                if 'Status=WAITING' in text:
                    return ('WAITING', None)
                elif 'Status=RUNNING' in text:
                    return ('RUNNING', None)
                elif 'Status=COMPLETE' in text:
                    return ('COMPLETE', text)
                else:
                    return ('FAILED', None)
            
            return ('FAILED', None)
        except Exception as e:
            print(f"检查状态失败: {e}")
            return ('FAILED', None)
    
    def run_blast(self, sequence: str, timeout: int = 300, poll_interval: int = 10) -> Optional[Dict]:
        """
        运行完整的BLAST流程
        
        参数:
            sequence: 查询序列
            timeout: 超时时间（秒）
            poll_interval: 轮询间隔（秒）
        
        返回:
            results: BLAST结果字典
        """
        request_id = self.submit_blast(sequence)
        
        if not request_id:
            return None
        
        start_time = time.time()
        
        while time.time() - start_time < timeout:
            status, result = self.check_status(request_id)
            
            if status == 'COMPLETE' and result:
                return self.parse_result(result)
            elif status == 'FAILED':
                return None
            
            time.sleep(poll_interval)
        
        return None
    
    def parse_result(self, xml_text: str) -> Dict:
        """
        解析BLAST XML结果
        
        返回格式:
        {
            'query_id': str,
            'query_length': int,
            'hits': [
                {
                    'accession': str,
                    'description': str,
                    'identity': float,
                    'length': int,
                    'e_value': float,
                    'score': float,
                    'strand': str,
                    'alignments': []
                }
            ]
        }
        """
        try:
            root = ET.fromstring(xml_text)
            
            results = {
                'query_id': '',
                'query_length': 0,
                'hits': []
            }
            
            query_info = root.find('.//Query-def')
            if query_info is not None:
                results['query_id'] = query_info.text
            
            query_len = root.find('.//Query-len')
            if query_len is not None:
                results['query_length'] = int(query_len.text)
            
            for hit in root.findall('.//Hit'):
                hit_info = {}
                
                accession = hit.find('.//Hit_accession')
                if accession is not None:
                    hit_info['accession'] = accession.text
                
                desc = hit.find('.//Hit_def')
                if desc is not None:
                    hit_info['description'] = desc.text
                
                hit_len = hit.find('.//Hit_len')
                if hit_len is not None:
                    hit_info['length'] = int(hit_len.text)
                
                for hsp in hit.findall('.//Hsp'):
                    identity = hsp.find('.//Hsp_identity')
                    if identity is not None:
                        hit_info['identity'] = int(identity.text)
                    
                    align_len = hsp.find('.//Hsp_align-len')
                    if align_len is not None:
                        align_length = int(align_len.text)
                        if 'identity' in hit_info and align_length > 0:
                            hit_info['identity_percent'] = (hit_info['identity'] / align_length) * 100
                    
                    e_value = hsp.find('.//Hsp_evalue')
                    if e_value is not None:
                        hit_info['e_value'] = float(e_value.text)
                    
                    score = hsp.find('.//Hsp_score')
                    if score is not None:
                        hit_info['score'] = float(score.text)
                    
                    strand = hsp.find('.//Hsp_strand')
                    if strand is not None:
                        hit_info['strand'] = strand.text
                
                if 'identity_percent' in hit_info:
                    results['hits'].append(hit_info)
            
            return results
        except Exception as e:
            print(f"解析结果失败: {e}")
            return None
    
    def get_off_target_risk(self, results: Dict, target_gene: str = "", threshold: float = 80.0) -> Dict:
        """
        分析脱靶风险
        
        参数:
            results: BLAST结果字典
            target_gene: 目标基因名称（用于排除自身）
            threshold: 相似度阈值（百分比）
        
        返回:
            risk分析结果
        """
        risk_result = {
            'has_off_target': False,
            'high_risk': [],
            'medium_risk': [],
            'low_risk': [],
            'total_hits': 0,
            'target_excluded': False
        }
        
        for hit in results.get('hits', []):
            identity = hit.get('identity_percent', 0)
            
            if identity >= threshold:
                is_target = False
                
                if target_gene:
                    desc = hit.get('description', '').lower()
                    if target_gene.lower() in desc:
                        is_target = True
                        risk_result['target_excluded'] = True
                
                if not is_target:
                    risk_result['has_off_target'] = True
                    
                    hit_info = {
                        'accession': hit.get('accession', ''),
                        'description': hit.get('description', ''),
                        'identity': identity,
                        'e_value': hit.get('e_value', 0),
                        'strand': hit.get('strand', '')
                    }
                    
                    if identity >= 90:
                        risk_result['high_risk'].append(hit_info)
                    elif identity >= 85:
                        risk_result['medium_risk'].append(hit_info)
                    else:
                        risk_result['low_risk'].append(hit_info)
        
        risk_result['total_hits'] = len(risk_result['high_risk']) + len(risk_result['medium_risk']) + len(risk_result['low_risk'])
        
        return risk_result

if __name__ == "__main__":
    blast = NCBIBLAST()
    
    test_sequence = "AUGAGGAGAUCUCCUCAGGCA"
    
    print("提交BLAST请求...")
    result = blast.run_blast(test_sequence)
    
    if result:
        print("\nBLAST结果:")
        print(f"查询序列: {result['query_id']}")
        print(f"序列长度: {result['query_length']}")
        print(f"命中数量: {len(result['hits'])}")
        
        for hit in result['hits'][:5]:
            print(f"\nAccession: {hit.get('accession')}")
            print(f"描述: {hit.get('description')}")
            print(f"相似度: {hit.get('identity_percent', 0):.2f}%")
            print(f"E-value: {hit.get('e_value')}")
    else:
        print("BLAST请求失败")